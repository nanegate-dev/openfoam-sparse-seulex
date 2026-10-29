/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
    Modified code Copyright (C) 2026 Nanegate
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::seulexFast

Group
    grpODESolvers

Description
    seulex (an extrapolation algorithm based on the linearly implicit Euler
    method with step size control and order selection) with the same
    algorithm, constants and step sequence, made faster without changing the
    answers beyond rounding:
      - the matrices (1/dx) I - J are factorised by a sparse LU on the
        sparsity pattern of the Jacobian instead of a dense LU
        (seulexSparseLU), whose forward substitution skips the leading zero
        block of the right-hand side;
      - with the switch autonomous (off by default), for ODE systems whose
        derivatives do not depend on the time argument, such as OpenFOAM's
        chemistry: the first right-hand side f(y0) of a step is evaluated
        once and shared by all extrapolation columns, and the right-hand
        side evaluated by the divergence test is reused for the next
        substep.

    The pattern is found from the Jacobian itself: at the first call from
    the Jacobian at the given state and at a state with every component
    shifted up (so entries that vanish only because a concentration is zero
    are included), and every later Jacobian is checked for entries outside
    it, which are added before the next factorisation. The sparse
    factorisation is used only when it needs fewer than half the
    multiply-adds of the dense one; a pivot that is too small without
    pivoting falls back to the dense LU with partial pivoting.

    A Jacobian that is a sparse matrix plus one outer product, J = S + u w^T,
    is split: OpenFOAM's chemistry in mass fractions adds rhoM*v_j*dYi/dt to
    every entry, which fills the whole matrix although S is as sparse as the
    concentration Jacobian. The split is found from the Jacobian (a robust
    median fit at the first call, then a cheap update on the entries outside
    the pattern), the sparse factorisation is of (1/dx) I - S, and the outer
    product is added back exactly in every solve by the Sherman-Morrison
    formula. When no split leaves the matrix sparse the solver works on J
    as before.

    Usage, in chemistryProperties:
    \verbatim
        ode
        {
            solver          seulexFast;
            absTol          1e-12;
            relTol          1e-1;
            autonomous      yes;    // chemistry: derivatives independent of t
        }
    \endverbatim

    Without autonomous (or with autonomous no) the solver makes every
    right-hand side evaluation seulex makes and is valid for any ODE system;
    only the sparse linear algebra is then used.

    Reference:
    \verbatim
        Hairer, E., Nørsett, S. P., & Wanner, G. (1996).
        Solving Ordinary Differential Equations II:
        Stiff and Differential-Algebraic Problems, second edition",
        Springer-Verlag, Berlin.
    \endverbatim

SourceFiles
    seulexFast.C

\*---------------------------------------------------------------------------*/

#ifndef seulexFast_H
#define seulexFast_H

#include "ODESolver.H"
#include "scalarMatrices.H"
#include "labelField.H"
#include "boolList.H"
#include "seulexSparseLU.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class seulexFast Declaration
\*---------------------------------------------------------------------------*/

class seulexFast
:
    public ODESolver
{
    // Private data

        // Static constants

            static const label kMaxx_ = 12;
            static const label iMaxx_ = kMaxx_ + 1;

            static const scalar
                stepFactor1_, stepFactor2_, stepFactor3_,
                stepFactor4_, stepFactor5_,
                kFactor1_, kFactor2_;

        //- The derivatives do not depend on the time argument: right-hand
        //  sides at the same state are shared (switch autonomous)
        const bool autonomous_;

        // Evaluated constants

            scalar jacRedo_;
            labelField nSeq_;
            scalarField cpu_;
            scalarSquareMatrix coeff_;

        // Temporary storage
        // held to avoid dynamic memory allocation between calls
        // and to transfer internal values between functions

            mutable scalar theta_;
            mutable label kTarg_;
            mutable scalarRectangularMatrix table_;

            mutable scalarField dfdx_;
            mutable scalarSquareMatrix dfdy_;
            mutable scalarSquareMatrix a_;
            mutable labelList pivotIndices_;

            // Fields space for "solve" function
            mutable scalarField dxOpt_, temp_;
            mutable scalarField y0_, ySequence_, scale_;

            // Fields used in "seul" function
            mutable scalarField dy_, yTemp_, dydx_;

            //- f(y0) of the current step attempt (autonomous systems)
            mutable scalarField f0_;

        // Sparse linear algebra

            //- Possible non-zeros of the Jacobian, n x n row-major
            mutable boolList pattern_;

            //- Ordering and symbolic factorisation for pattern_
            mutable seulexSparseLU slu_;

            //- The sparse factorisation is cheaper for this pattern
            mutable bool sparse_;

            //- The current factorisation is the dense one
            mutable bool denseNow_;

            //- Pattern rebuilds; past a limit (the system keeps changing
            //  size, as with TDAC) the solver stays dense
            mutable label nBuild_;

        // Split J = S + u w^T

            //- The Jacobian has a useful split (decided at the first call)
            mutable bool rankOne_;

            //- dfdy_ currently holds S rather than J
            mutable bool split_;

            //- The outer product's vectors
            mutable scalarField u_, w_;

            //- (1/dx I - S)^-1 u of the current sparse factorisation and
            //  1 - w.z
            mutable scalarField z_;
            mutable scalar smDenom_;

            //- Columns with few pattern entries, used to update u
            mutable labelList uCols_;

            //- Work space of the medians
            mutable scalarList ratios_;


    // Private Member Functions

        //- Pattern rebuilds allowed: entries below rounding against the
        //  outer product at the first state appear later, so a split
        //  Jacobian is allowed more
        label maxBuild() const
        {
            return rankOne_ ? 200 : 10;
        }

        //- Median of the first n entries of ratios_ (reordered)
        scalar median(const label n) const;

        //- Robust fit of J's outer product into u, w: for the columns with
        //  the most non-zeros as reference, w from median column ratios and
        //  u from median row ratios; returns the number of entries the best
        //  fit leaves
        label fitRankOne
        (
            const scalarSquareMatrix& J,
            scalarField& u,
            scalarField& w
        ) const;

        //- Update u, w from the entries of J outside the pattern, starting
        //  from the current w
        void updateRankOne(const scalarSquareMatrix& J) const;

        //- Replace J by J - u w^T, zeroing what is left by rounding only;
        //  returns the number of non-zeros left outside the pattern (all
        //  non-zeros when the pattern is not yet set)
        label subtractRankOne(scalarSquareMatrix& J) const;

        //- Split dfdy_ into S and u w^T when the Jacobian has a useful split
        void splitJacobian() const;

        //- Add the Jacobian's non-zeros to the pattern; rebuild the
        //  symbolic factorisation when it grew
        void updatePattern
        (
            const scalar x,
            const scalarField& y,
            const label li
        ) const;

        //- Factorise (1/dx) I - J
        void factorise(const scalar dx) const;

        //- Solve with the current factorisation
        void linSolve(scalarField& b) const;

        //- Computes the j-th line of the extrapolation table
        bool seul
        (
            const scalar x0,
            const scalarField& y0,
            const label li,
            const scalar dxTot,
            const label k,
            scalarField& y,
            const scalarField& scale
        ) const;

        //- Polynomial extrapolation
        void extrapolate
        (
            const label k,
            scalarRectangularMatrix& table,
            scalarField& y
        ) const;


public:

    //- Runtime type information
    TypeName("seulexFast");


    // Constructors

        //- Construct from ODESystem
        seulexFast(const ODESystem& ode, const dictionary& dict);


    //- Destructor
    virtual ~seulexFast() = default;


    // Member Functions

        //- Resize the ODE solver
        virtual bool resize();

        //- Solve the ODE system and the update the state
        virtual void solve
        (
            scalar& x,
            scalarField& y,
            const label li,
            stepState& step
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
