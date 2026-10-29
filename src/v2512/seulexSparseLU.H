/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Nanegate
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
    Foam::seulexSparseLU

Description
    LU factorisation of A = r*I - J on a fixed sparsity pattern of the
    Jacobian J, for the linear systems of seulexFast.

    The ordering (minimum degree on the symmetrised pattern) and the symbolic
    factorisation are computed once per pattern; a numeric factorisation is
    then a flat list of multiply-adds over the non-zeros and the fill-in
    only. There is no pivoting: factor() returns false when a pivot is tiny
    against the diagonal it started from (or not finite), and the caller then
    uses the dense LU with partial pivoting for that matrix.

    For chemistry the Jacobian is very sparse once a mechanism has more than
    a few dozen species (4 % of the entries are non-zero for a 544-species
    n-heptane mechanism), so the factorisation costs a small fraction of the
    dense one.

SourceFiles
    seulexSparseLU.C
    seulexSparseLUI.H

\*---------------------------------------------------------------------------*/

#ifndef seulexSparseLU_H
#define seulexSparseLU_H

#include "scalarMatrices.H"
#include "boolList.H"
#include "labelList.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class seulexSparseLU Declaration
\*---------------------------------------------------------------------------*/

class seulexSparseLU
{
    // Private Data

        //- Number of unknowns
        label n_;

        //- Elimination order: row i of the factors is row perm_[i] of A
        labelList perm_;

        //- Inverse of perm_
        labelList iperm_;

        //- Row starts of the stored L+U entries (in the elimination order)
        labelList rowStart_;

        //- Column of each stored entry (in the elimination order)
        labelList col_;

        //- Position of the diagonal entry of each row
        labelList diag_;

        //- Position of the first upper-triangular entry of each row
        labelList upper_;

        //- Position of each stored entry in the dense J; -1 for fill-in
        labelList src_;

        //- Elimination program: the eliminated entry, its pivot and the
        //  number of updates of each elimination ...
        labelList elimP_;
        labelList elimPiv_;
        labelList elimCnt_;

        //- ... and the destinations and sources of all the updates, in order
        labelList updDst_;
        labelList updSrc_;

        //- Values of the factors
        scalarList val_;

        //- Pivot scale of each row: its diagonal before elimination
        scalarList scale_;

        //- Work space of the solves
        scalarList work_;

        //- Multiply-adds of one numeric factorisation
        scalar nUpdates_;


public:

    // Constructors

        //- Construct empty; build() sets the pattern
        seulexSparseLU();


    // Member Functions

        //- Compute the ordering and the symbolic factorisation for the
        //  possible non-zeros of J (n x n flags, row-major)
        void build(const label n, const boolList& pattern);

        //- Number of unknowns
        label n() const
        {
            return n_;
        }

        //- Multiply-adds of one numeric factorisation
        scalar nUpdates() const
        {
            return nUpdates_;
        }

        //- Take entry d (row-major in the dense J) into the factorisation
        //  when the factors already store its position (fill-in); false
        //  when the pattern has to be rebuilt for it
        bool addEntry(const label d);

        //- Factorise r*I - J; false if a pivot is too small without
        //  pivoting (the caller then factorises densely)
        inline bool factor(const scalarSquareMatrix& J, const scalar r);

        //- Solve A x = b in place with the current factors
        inline void solve(scalarField& b);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "seulexSparseLUI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
