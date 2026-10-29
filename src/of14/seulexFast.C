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

\*---------------------------------------------------------------------------*/

#include "seulexFast.H"
#include "addToRunTimeSelectionTable.H"
#include "Switch.H"
#include "DynamicList.H"

#include <algorithm>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(seulexFast, 0);

    addToRunTimeSelectionTable(ODESolver, seulexFast, dictionary);

    const scalar
        seulexFast::stepFactor1_ = 0.6,
        seulexFast::stepFactor2_ = 0.93,
        seulexFast::stepFactor3_ = 0.1,
        seulexFast::stepFactor4_ = 4,
        seulexFast::stepFactor5_ = 0.5,
        seulexFast::kFactor1_ = 0.7,
        seulexFast::kFactor2_ = 0.9;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::seulexFast::seulexFast(const ODESystem& ode, const dictionary& dict)
:
    ODESolver(ode, dict),
    autonomous_(dict.lookupOrDefault<Switch>("autonomous", false)),
    jacRedo_(min(1e-4, min(relTol_))),
    nSeq_(iMaxx_),
    cpu_(iMaxx_),
    coeff_(iMaxx_, iMaxx_),
    theta_(2*jacRedo_),
    table_(kMaxx_, n_),
    dfdx_(n_),
    dfdy_(n_),
    a_(n_),
    pivotIndices_(n_),
    dxOpt_(iMaxx_),
    temp_(iMaxx_),
    y0_(n_),
    ySequence_(n_),
    scale_(n_),
    dy_(n_),
    yTemp_(n_),
    dydx_(n_),
    f0_(n_),
    sparse_(false),
    denseNow_(true),
    nBuild_(0),
    rankOne_(false),
    split_(false),
    u_(n_, 0),
    w_(n_, 0),
    z_(n_, 0),
    smDenom_(1),
    ratios_(n_)
{
    // The CPU time factors for the major parts of the algorithm
    const scalar cpuFunc = 1, cpuJac = 5, cpuLU = 1, cpuSolve = 1;

    nSeq_[0] = 2;
    nSeq_[1] = 3;

    for (int i=2; i<iMaxx_; i++)
    {
        nSeq_[i] = 2*nSeq_[i-2];
    }
    cpu_[0] = cpuJac + cpuLU + nSeq_[0]*(cpuFunc + cpuSolve);

    for (int k=0; k<kMaxx_; k++)
    {
        cpu_[k+1] = cpu_[k] + (nSeq_[k+1]-1)*(cpuFunc + cpuSolve) + cpuLU;
    }

    // Set the extrapolation coefficients array
    for (int k=0; k<iMaxx_; k++)
    {
        for (int l=0; l<k; l++)
        {
            scalar ratio = scalar(nSeq_[k])/nSeq_[l];
            coeff_(k, l) = 1/(ratio - 1);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::scalar Foam::seulexFast::median(const label n) const
{
    scalar* b = ratios_.begin();
    std::nth_element(b, b + n/2, b + n);
    return b[n/2];
}


Foam::label Foam::seulexFast::fitRankOne
(
    const scalarSquareMatrix& J,
    scalarField& u,
    scalarField& w
) const
{
    const label n = n_;

    // Candidate reference columns: the ones with the most non-zeros
    labelList nnz(n, 0);
    for (label i = 0; i < n; ++i)
    {
        for (label j = 0; j < n; ++j)
        {
            if (J(i, j) != 0)
            {
                ++nnz[j];
            }
        }
    }

    const label nCand = min(label(8), n);
    labelList cand(nCand, -1);
    for (label c = 0; c < nCand; ++c)
    {
        label best = -1;
        for (label j = 0; j < n; ++j)
        {
            if (nnz[j] >= 0 && (best < 0 || nnz[j] > nnz[best]))
            {
                best = j;
            }
        }
        cand[c] = best;
        nnz[best] = -1;
    }

    scalarField uk(n), wk(n);
    label bestLeft = labelMax;

    forAll(cand, c)
    {
        const label k = cand[c];

        // w: the median ratio of each column to column k
        for (label j = 0; j < n; ++j)
        {
            label m = 0;
            for (label i = 0; i < n; ++i)
            {
                if (mag(J(i, k)) > rootVSmall)
                {
                    ratios_[m++] = J(i, j)/J(i, k);
                }
            }
            wk[j] = m ? median(m) : 0;
        }

        // u: the median ratio of each row to w
        for (label i = 0; i < n; ++i)
        {
            label m = 0;
            for (label j = 0; j < n; ++j)
            {
                if (mag(wk[j]) > rootVSmall)
                {
                    ratios_[m++] = J(i, j)/wk[j];
                }
            }
            uk[i] = m ? median(m) : 0;
        }

        label left = 0;
        for (label i = 0; i < n; ++i)
        {
            for (label j = 0; j < n; ++j)
            {
                const scalar uw = uk[i]*wk[j];
                if (mag(J(i, j) - uw) > 1e-12*max(mag(J(i, j)), mag(uw)))
                {
                    ++left;
                }
            }
        }

        if (left < bestLeft)
        {
            bestLeft = left;
            u = uk;
            w = wk;
        }
    }

    return bestLeft;
}


void Foam::seulexFast::updateRankOne(const scalarSquareMatrix& J) const
{
    const label n = n_;
    const label nRows = uCols_.size();

    if (nRows == 0)
    {
        u_ = 0;
        w_ = 0;
        return;
    }

    // u from the entries outside the pattern of the columns uCols_; there
    // J(i, k) = u_i w_k to rounding
    auto fitU = [&]()
    {
        for (label i = 0; i < n; ++i)
        {
            label m = 0;
            forAll(uCols_, c)
            {
                const label k = uCols_[c];
                if (!pattern_[i*n + k] && mag(w_[k]) > rootVSmall)
                {
                    ratios_[m++] = J(i, k)/w_[k];
                }
            }
            u_[i] = m ? median(m) : 0;
        }
    };

    fitU();

    // w from the entries outside the pattern of the rows with the largest
    // |u|
    labelList rows(nRows, -1);
    label nr = 0;
    {
        for (label i = 0; i < n; ++i)
        {
            ratios_[i] = mag(u_[i]);
        }
        const label kth = min(n - 1, max(n - nRows, label(0)));
        scalar* b = ratios_.begin();
        std::nth_element(b, b + kth, b + n);
        const scalar threshold = b[kth];
        for (label i = 0; i < n && nr < nRows; ++i)
        {
            if (mag(u_[i]) > rootVSmall && mag(u_[i]) >= threshold)
            {
                rows[nr++] = i;
            }
        }
    }

    for (label j = 0; j < n; ++j)
    {
        label m = 0;
        for (label r = 0; r < nr; ++r)
        {
            const label i = rows[r];
            if (!pattern_[i*n + j])
            {
                ratios_[m++] = J(i, j)/u_[i];
            }
        }
        w_[j] = m ? median(m) : 0;
    }

    fitU();
}


Foam::label Foam::seulexFast::subtractRankOne(scalarSquareMatrix& J) const
{
    const label n = n_;
    const bool havePattern = pattern_.size() == n*n;
    label outside = 0;

    for (label i = 0; i < n; ++i)
    {
        for (label j = 0; j < n; ++j)
        {
            scalar& Jij = J(i, j);
            const scalar uw = u_[i]*w_[j];
            if (uw != 0)
            {
                const scalar s = Jij - uw;
                Jij = mag(s) > 1e-12*max(mag(Jij), mag(uw)) ? s : 0;
            }
            if (Jij != 0 && (!havePattern || !pattern_[i*n + j]))
            {
                ++outside;
            }
        }
    }

    return outside;
}


void Foam::seulexFast::splitJacobian() const
{
    split_ = false;

    if (!rankOne_)
    {
        return;
    }

    // J is kept in a_ (free until the next factorisation) in case no split
    // fits
    a_ = dfdy_;

    updateRankOne(dfdy_);
    label outside = subtractRankOne(dfdy_);

    if (outside > n_)
    {
        // The update did not fit (the outer product turned or the pattern
        // misses entries): fit again from scratch
        dfdy_ = a_;
        fitRankOne(dfdy_, u_, w_);
        outside = subtractRankOne(dfdy_);

        if (debug)
        {
            Info<< type() << ": split fitted again, " << outside
                << " entries outside the pattern" << endl;
        }
    }

    if (outside > n_)
    {
        // No useful split any more: work on J
        dfdy_ = a_;
        rankOne_ = false;
        sparse_ = false;
        u_ = 0;
        w_ = 0;

        if (debug)
        {
            Info<< type() << ": split switched off" << endl;
        }
        return;
    }

    split_ = true;
}


void Foam::seulexFast::updatePattern
(
    const scalar x,
    const scalarField& y,
    const label li
) const
{
    const label n = n_;

    if (nBuild_ > maxBuild() || (!sparse_ && pattern_.size() == n*n))
    {
        sparse_ = false;
        split_ = false;
        return;
    }

    bool grown = false;

    if (pattern_.size() != n*n)
    {
        if (nBuild_ >= maxBuild())
        {
            ++nBuild_;
            sparse_ = false;
            split_ = false;
            return;
        }

        // Entries that vanish only at this state (zero concentrations):
        // the Jacobian at a state with every component shifted up
        scalarField yProbe(y);
        forAll(yProbe, i)
        {
            const scalar yi = y[i];
            yProbe[i] =
                yi >= -1e-6
              ? mag(yi) + 1e-6*(1 + mag(yi))
              : yi - 1e-6*(1 + mag(yi));
        }
        scalarField dfdxProbe(n);
        scalarSquareMatrix dfdyProbe(n);
        odes_.jacobian(x, yProbe, li, dfdxProbe, dfdyProbe);

        // Split off an outer product when that leaves the probe Jacobian
        // at most a quarter as full
        label nnz = 0;
        const scalar* p = dfdyProbe.v();
        for (label q = 0; q < n*n; ++q)
        {
            if (p[q] != 0)
            {
                ++nnz;
            }
        }
        rankOne_ = fitRankOne(dfdyProbe, u_, w_) < nnz/4;
        if (rankOne_)
        {
            subtractRankOne(dfdyProbe);
        }
        else
        {
            u_ = 0;
            w_ = 0;
        }

        pattern_.resize(n*n);
        pattern_ = false;
        for (label i = 0; i < n; ++i)
        {
            pattern_[i*n + i] = true;
        }
        forAll(pattern_, q)
        {
            if (p[q] != 0)
            {
                pattern_[q] = true;
            }
        }

        if (rankOne_)
        {
            // The columns with the fewest pattern entries update u
            labelList nPat(n, 0);
            for (label i = 0; i < n; ++i)
            {
                for (label j = 0; j < n; ++j)
                {
                    if (pattern_[i*n + j])
                    {
                        ++nPat[j];
                    }
                }
            }
            DynamicList<label> cols;
            for (label c = 0; c < min(label(16), n); ++c)
            {
                label best = -1;
                for (label j = 0; j < n; ++j)
                {
                    if
                    (
                        w_[j] != 0 && nPat[j] >= 0
                     && (best < 0 || nPat[j] < nPat[best])
                    )
                    {
                        best = j;
                    }
                }
                if (best < 0)
                {
                    break;
                }
                cols.append(best);
                nPat[best] = -1;
            }
            uCols_.transfer(cols);
        }

        grown = true;
    }

    splitJacobian();

    if (nBuild_ > maxBuild())
    {
        if (split_)
        {
            dfdy_ = a_;
            split_ = false;
        }
        rankOne_ = false;
        sparse_ = false;
        return;
    }

    // New non-zeros: taken into the current factors where those already
    // store the position (fill-in), otherwise the pattern is rebuilt
    const scalar* J = dfdy_.v();
    label nTaken = 0;
    forAll(pattern_, q)
    {
        if (J[q] != 0 && !pattern_[q])
        {
            pattern_[q] = true;
            if (!grown && slu_.addEntry(q))
            {
                ++nTaken;
            }
            else
            {
                grown = true;
            }
        }
    }

    if (debug && nTaken && !grown)
    {
        Info<< type() << ": " << nTaken
            << " new entries taken into the fill-in" << endl;
    }

    if (grown)
    {
        if (++nBuild_ > maxBuild())
        {
            if (split_)
            {
                dfdy_ = a_;
                split_ = false;
            }
            rankOne_ = false;
            sparse_ = false;
            return;
        }

        slu_.build(n, pattern_);
        sparse_ = slu_.nUpdates() < 0.5*scalar(n)*n*n/3;
        if (!sparse_ && split_)
        {
            dfdy_ = a_;
            split_ = false;
            rankOne_ = false;
        }

        if (debug)
        {
            Info<< type() << ": pattern built (" << nBuild_ << "), n " << n
                << ", split " << (rankOne_ ? "yes" : "no")
                << ", sparse factorisation " << slu_.nUpdates()
                << " multiply-adds against " << scalar(n)*n*n/3
                << " dense, using " << (sparse_ ? "sparse" : "dense") << endl;
        }
    }
}


void Foam::seulexFast::factorise(const scalar dx) const
{
    denseNow_ = !(sparse_ && slu_.factor(dfdy_, 1/dx));

    if (!denseNow_ && split_)
    {
        // (1/dx I - S - u w^T)^-1 b = B^-1 b + z (w.B^-1 b)/(1 - w.z) with
        // B = 1/dx I - S and z = B^-1 u (Sherman-Morrison)
        z_ = u_;
        slu_.solve(z_);
        scalar wz = 0;
        forAll(w_, i)
        {
            wz += w_[i]*z_[i];
        }
        smDenom_ = 1 - wz;
        denseNow_ = !(mag(smDenom_) > 1e-8*(1 + mag(wz)));
    }

    if (denseNow_)
    {
        for (label i=0; i<n_; i++)
        {
            for (label j=0; j<n_; j++)
            {
                a_(i, j) = -dfdy_(i, j);
            }

            if (split_)
            {
                for (label j=0; j<n_; j++)
                {
                    a_(i, j) -= u_[i]*w_[j];
                }
            }

            a_(i, i) += 1/dx;
        }

        LUDecompose(a_, pivotIndices_);
    }
}


void Foam::seulexFast::linSolve(scalarField& b) const
{
    if (denseNow_)
    {
        LUBacksubstitute(a_, pivotIndices_, b);
    }
    else
    {
        slu_.solve(b);

        if (split_)
        {
            scalar wb = 0;
            forAll(w_, i)
            {
                wb += w_[i]*b[i];
            }
            const scalar c = wb/smDenom_;
            forAll(b, i)
            {
                b[i] += c*z_[i];
            }
        }
    }
}


bool Foam::seulexFast::seul
(
    const scalar x0,
    const scalarField& y0,
    const label li,
    const scalar dxTot,
    const label k,
    scalarField& y,
    const scalarField& scale
) const
{
    label nSteps = nSeq_[k];
    scalar dx = dxTot/nSteps;

    factorise(dx);

    scalar xnew = x0 + dx;

    if (autonomous_)
    {
        // The derivatives do not depend on the time argument, so every
        // column's first right-hand side f(x0 + dx_k, y0) is the same
        // vector: evaluated with the first column of each step attempt
        if (k == 0)
        {
            odes_.derivatives(xnew, y0, li, f0_);
        }
        dy_ = f0_;
    }
    else
    {
        odes_.derivatives(xnew, y0, li, dy_);
    }
    linSolve(dy_);

    yTemp_ = y0;

    // dydx_ holds f(yTemp_) after the divergence test (autonomous systems)
    bool haveF = false;

    for (label nn=1; nn<nSteps; nn++)
    {
        yTemp_ += dy_;
        xnew += dx;

        if (nn == 1 && k<=1)
        {
            scalar dy1 = 0;
            for (label i=0; i<n_; i++)
            {
                dy1 += sqr(dy_[i]/scale[i]);
            }
            dy1 = sqrt(dy1);

            odes_.derivatives(x0 + dx, yTemp_, li, dydx_);
            haveF = autonomous_;
            for (label i=0; i<n_; i++)
            {
                dy_[i] = dydx_[i] - dy_[i]/dx;
            }

            linSolve(dy_);

            const scalar denom = max(1, dy1);
            scalar dy2 = 0;
            for (label i=0; i<n_; i++)
            {
                // Test of dy_[i] to avoid overflow
                if (mag(dy_[i]) > scale[i]*denom)
                {
                    theta_ = 1;
                    return false;
                }

                dy2 += sqr(dy_[i]/scale[i]);
            }
            dy2 = sqrt(dy2);
            theta_ = dy2/denom;

            if (theta_ > 1)
            {
                return false;
            }
        }

        if (haveF)
        {
            // The divergence test has just evaluated f at this same state
            dy_ = dydx_;
            haveF = false;
        }
        else
        {
            odes_.derivatives(xnew, yTemp_, li, dy_);
        }
        linSolve(dy_);
    }

    for (label i=0; i<n_; i++)
    {
        y[i] = yTemp_[i] + dy_[i];
    }

    return true;
}


void Foam::seulexFast::extrapolate
(
    const label k,
    scalarRectangularMatrix& table,
    scalarField& y
) const
{
    for (int j=k-1; j>0; j--)
    {
        for (label i=0; i<n_; i++)
        {
            table[j-1][i] =
                table(j, i) + coeff_(k, j)*(table(j, i) - table[j-1][i]);
        }
    }

    for (int i=0; i<n_; i++)
    {
        y[i] = table(0, i) + coeff_(k, 0)*(table(0, i) - y[i]);
    }
}


bool Foam::seulexFast::resize()
{
    if (ODESolver::resize())
    {
        table_.shallowResize(kMaxx_, n_);
        resizeField(dfdx_);
        resizeMatrix(dfdy_);
        resizeMatrix(a_);
        resizeField(pivotIndices_);
        resizeField(y0_);
        resizeField(ySequence_);
        resizeField(scale_);
        resizeField(dy_);
        resizeField(yTemp_);
        resizeField(dydx_);
        resizeField(f0_);
        resizeField(u_);
        resizeField(w_);
        resizeField(z_);
        ratios_.setSize(n_);

        pattern_.clear();
        sparse_ = false;
        denseNow_ = true;
        rankOne_ = false;
        split_ = false;
        uCols_.clear();

        return true;
    }

    return false;
}


void Foam::seulexFast::solve
(
    scalar& x,
    scalarField& y,
    const label li,
    stepState& step
) const
{
    temp_[0] = great;
    scalar dx = step.dxTry;
    y0_ = y;
    dxOpt_[0] = mag(0.1*dx);

    if (step.first || step.prevReject)
    {
        theta_ = 2*jacRedo_;
    }

    if (step.first)
    {
        // NOTE: the first element of relTol_ and absTol_ are used here.
        scalar logTol = -log10(relTol_[0] + absTol_[0])*0.6 + 0.5;
        kTarg_ = max(1, min(kMaxx_ - 1, int(logTol)));
    }

    forAll(scale_, i)
    {
        scale_[i] = absTol_[i] + relTol_[i]*mag(y[i]);
    }

    bool jacUpdated = false;

    if (theta_ > jacRedo_)
    {
        odes_.jacobian(x, y, li, dfdx_, dfdy_);
        updatePattern(x, y, li);
        jacUpdated = true;
    }

    int k;
    scalar dxNew = mag(dx);
    bool firstk = true;

    while (firstk || step.reject)
    {
        dx = step.forward ? dxNew : -dxNew;
        firstk = false;
        step.reject = false;

        if (mag(dx) <= mag(x)*sqr(small))
        {
             WarningInFunction
                    << "step size underflow :"  << dx << endl;
        }

        scalar errOld = 0;

        for (k=0; k<=kTarg_+1; k++)
        {
            bool success = seul(x, y0_, li, dx, k, ySequence_, scale_);

            if (!success)
            {
                step.reject = true;
                dxNew = mag(dx)*stepFactor5_;
                break;
            }

            if (k == 0)
            {
                 y = ySequence_;
            }
            else
            {
                forAll(ySequence_, i)
                {
                    table_[k-1][i] = ySequence_[i];
                }
            }

            if (k != 0)
            {
                extrapolate(k, table_, y);
                scalar err = 0;
                forAll(scale_, i)
                {
                    scale_[i] = absTol_[i] + relTol_[i]*mag(y0_[i]);
                    err += sqr((y[i] - table_(0, i))/scale_[i]);
                }
                err = sqrt(err/n_);
                if (err > 1/small || (k > 1 && err >= errOld))
                {
                    step.reject = true;
                    dxNew = mag(dx)*stepFactor5_;
                    break;
                }
                errOld = max(4*err, 1);
                scalar expo = 1.0/(k + 1);
                scalar facmin = pow(stepFactor3_, expo);
                scalar fac;
                if (err == 0)
                {
                    fac = 1/facmin;
                }
                else
                {
                    fac = stepFactor2_/pow(err/stepFactor1_, expo);
                    fac = max(facmin/stepFactor4_, min(1/facmin, fac));
                }
                dxOpt_[k] = mag(dx*fac);
                temp_[k] = cpu_[k]/dxOpt_[k];

                if ((step.first || step.last) && err <= 1)
                {
                    break;
                }

                if
                (
                    k == kTarg_ - 1
                 && !step.prevReject
                 && !step.first && !step.last
                )
                {
                    if (err <= 1)
                    {
                        break;
                    }
                    else if (err > nSeq_[kTarg_]*nSeq_[kTarg_ + 1]*4)
                    {
                        step.reject = true;
                        kTarg_ = k;
                        if (kTarg_>1 && temp_[k-1] < kFactor1_*temp_[k])
                        {
                            kTarg_--;
                        }
                        dxNew = dxOpt_[kTarg_];
                        break;
                    }
                }

                if (k == kTarg_)
                {
                    if (err <= 1)
                    {
                        break;
                    }
                    else if (err > nSeq_[k + 1]*2)
                    {
                        step.reject = true;
                        if (kTarg_>1 && temp_[k-1] < kFactor1_*temp_[k])
                        {
                            kTarg_--;
                        }
                        dxNew = dxOpt_[kTarg_];
                        break;
                    }
                }

                if (k == kTarg_+1)
                {
                    if (err > 1)
                    {
                        step.reject = true;
                        if
                        (
                            kTarg_ > 1
                         && temp_[kTarg_-1] < kFactor1_*temp_[kTarg_]
                        )
                        {
                            kTarg_--;
                        }
                        dxNew = dxOpt_[kTarg_];
                    }
                    break;
                }
            }
        }
        if (step.reject)
        {
            step.prevReject = true;
            if (!jacUpdated)
            {
                theta_ = 2*jacRedo_;

                if (theta_ > jacRedo_ && !jacUpdated)
                {
                    odes_.jacobian(x, y, li, dfdx_, dfdy_);
                    updatePattern(x, y, li);
                    jacUpdated = true;
                }
            }
        }
    }

    jacUpdated = false;

    step.dxDid = dx;
    x += dx;

    label kopt;
    if (k == 1)
    {
        kopt = 2;
    }
    else if (k <= kTarg_)
    {
        kopt=k;
        if (temp_[k-1] < kFactor1_*temp_[k])
        {
            kopt = k - 1;
        }
        else if (temp_[k] < kFactor2_*temp_[k - 1])
        {
            kopt = min(k + 1, kMaxx_ - 1);
        }
    }
    else
    {
        kopt = k - 1;
        if (k > 2 && temp_[k-2] < kFactor1_*temp_[k - 1])
        {
            kopt = k - 2;
        }
        if (temp_[k] < kFactor2_*temp_[kopt])
        {
            kopt = min(k, kMaxx_ - 1);
        }
    }

    if (step.prevReject)
    {
        kTarg_ = min(kopt, k);
        dxNew = min(mag(dx), dxOpt_[kTarg_]);
        step.prevReject = false;
    }
    else
    {
        if (kopt <= k)
        {
            dxNew = dxOpt_[kopt];
        }
        else
        {
            if (k < kTarg_ && temp_[k] < kFactor2_*temp_[k - 1])
            {
                dxNew = dxOpt_[k]*cpu_[kopt + 1]/cpu_[k];
            }
            else
            {
                dxNew = dxOpt_[k]*cpu_[kopt]/cpu_[k];
            }
        }
        kTarg_ = kopt;
    }

    step.dxTry = step.forward ? dxNew : -dxNew;
}


// ************************************************************************* //
