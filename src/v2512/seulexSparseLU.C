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

\*---------------------------------------------------------------------------*/

#include "seulexSparseLU.H"
#include "DynamicList.H"

#include <algorithm>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::seulexSparseLU::seulexSparseLU()
:
    n_(0),
    nUpdates_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::seulexSparseLU::build(const label n, const boolList& pattern)
{
    n_ = n;
    const scalar maxUpdates = 0.5*scalar(n)*n*n/3;

    label nnz = 0;
    forAll(pattern, q)
    {
        if (pattern[q])
        {
            ++nnz;
        }
    }
    if (nnz > 0.5*scalar(n)*n)
    {
        nUpdates_ = maxUpdates;
        rowStart_.clear();
        col_.clear();
        diag_.clear();
        upper_.clear();
        src_.clear();
        elimP_.clear();
        elimPiv_.clear();
        elimCnt_.clear();
        updDst_.clear();
        updSrc_.clear();
        val_.clear();
        return;
    }

    // Symmetrised adjacency of the pattern, without the diagonal
    List<boolList> adj(n, boolList(n, false));
    for (label i = 0; i < n; ++i)
    {
        for (label j = 0; j < n; ++j)
        {
            if (i != j && pattern[i*n + j])
            {
                adj[i][j] = true;
                adj[j][i] = true;
            }
        }
    }

    // Minimum degree ordering (greedy, ties to the lowest index)
    labelList degree(n, Zero);
    for (label i = 0; i < n; ++i)
    {
        for (label j = 0; j < n; ++j)
        {
            if (adj[i][j])
            {
                ++degree[i];
            }
        }
    }

    boolList done(n, false);
    perm_.resize(n);
    DynamicList<label> nb(n);

    for (label step = 0; step < n; ++step)
    {
        label best = -1;
        for (label i = 0; i < n; ++i)
        {
            if (!done[i] && (best < 0 || degree[i] < degree[best]))
            {
                best = i;
            }
        }
        perm_[step] = best;
        done[best] = true;

        // If the minimum degree equals n - 1 - step, the remaining vertices
        // already form a complete clique and need no further edge updates
        if (degree[best] == n - 1 - step)
        {
            for (label i = 0; i < n; ++i)
            {
                if (!done[i])
                {
                    perm_[++step] = i;
                    done[i] = true;
                }
            }
            break;
        }

        nb.clear();
        for (label j = 0; j < n; ++j)
        {
            if (!done[j] && adj[best][j])
            {
                nb.append(j);
            }
        }

        // Eliminating best connects all of its remaining neighbours
        for (const label a : nb)
        {
            adj[a][best] = false;
            --degree[a];
            if (degree[a] == n - 2 - step)
            {
                continue;
            }
            for (const label b : nb)
            {
                if (a != b && !adj[a][b])
                {
                    adj[a][b] = true;
                    ++degree[a];
                }
            }
        }
    }

    iperm_.resize(n);
    forAll(perm_, i)
    {
        iperm_[perm_[i]] = i;
    }

    // Symbolic factorisation, row by row in the new order
    nUpdates_ = 0;

    List<DynamicList<label>> rows(n);
    labelList upperStart(n, Zero);
    boolList mark(n);

    for (label i = 0; i < n; ++i)
    {
        mark = false;
        mark[i] = true;
        const label oi = perm_[i];
        for (label j = 0; j < n; ++j)
        {
            if (pattern[oi*n + perm_[j]])
            {
                mark[j] = true;
            }
        }
        for (label k = 0; k < i; ++k)
        {
            if (mark[k])
            {
                nUpdates_ += rows[k].size() - upperStart[k];
                for (label q = upperStart[k]; q < rows[k].size(); ++q)
                {
                    mark[rows[k][q]] = true;
                }
            }
        }
        if (nUpdates_ >= maxUpdates)
        {
            rowStart_.clear();
            col_.clear();
            diag_.clear();
            upper_.clear();
            src_.clear();
            elimP_.clear();
            elimPiv_.clear();
            elimCnt_.clear();
            updDst_.clear();
            updSrc_.clear();
            val_.clear();
            return;
        }
        for (label j = 0; j < n; ++j)
        {
            if (mark[j])
            {
                rows[i].append(j);
            }
        }
        upperStart[i] = label
        (
            std::lower_bound(rows[i].begin(), rows[i].end(), i + 1)
          - rows[i].begin()
        );
    }

    // Flat storage of L+U
    rowStart_.resize(n + 1);
    rowStart_[0] = 0;
    for (label i = 0; i < n; ++i)
    {
        rowStart_[i + 1] = rowStart_[i] + rows[i].size();
    }

    col_.resize(rowStart_[n]);
    diag_.resize(n);
    upper_.resize(n);
    for (label i = 0; i < n; ++i)
    {
        forAll(rows[i], q)
        {
            col_[rowStart_[i] + q] = rows[i][q];
        }
        diag_[i] = rowStart_[i] + label
        (
            std::lower_bound(rows[i].begin(), rows[i].end(), i)
          - rows[i].begin()
        );
        upper_[i] = rowStart_[i] + upperStart[i];
    }

    val_.resize(col_.size());
    val_ = 0;

    // Where each stored entry comes from in the dense J (-1: fill-in)
    src_.resize(col_.size());
    src_ = -1;
    for (label i = 0; i < n; ++i)
    {
        for (label p = rowStart_[i]; p < rowStart_[i + 1]; ++p)
        {
            const label d = perm_[i]*n + perm_[col_[p]];
            if (pattern[d])
            {
                src_[p] = d;
            }
        }
    }

    // Elimination program: for row i and each k < i in it,
    //   l = a_ik/u_kk; a_ik = l; a_ij -= l*u_kj for j in U(k)
    DynamicList<label> elimP, elimPiv, elimCnt, updDst, updSrc;
    labelList where(n, -1);
    nUpdates_ = 0;

    for (label i = 0; i < n; ++i)
    {
        for (label p = rowStart_[i]; p < rowStart_[i + 1]; ++p)
        {
            where[col_[p]] = p;
        }
        for (label p = rowStart_[i]; p < diag_[i]; ++p)
        {
            const label k = col_[p];
            const label cnt = rowStart_[k + 1] - upper_[k];
            elimP.append(p);
            elimPiv.append(diag_[k]);
            elimCnt.append(cnt);
            for (label q = upper_[k]; q < rowStart_[k + 1]; ++q)
            {
                updDst.append(where[col_[q]]);
                updSrc.append(q);
            }
            nUpdates_ += cnt;
        }
        for (label p = rowStart_[i]; p < rowStart_[i + 1]; ++p)
        {
            where[col_[p]] = -1;
        }
    }

    elimP_.transfer(elimP);
    elimPiv_.transfer(elimPiv);
    elimCnt_.transfer(elimCnt);
    updDst_.transfer(updDst);
    updSrc_.transfer(updSrc);

    scale_.resize(n);
    work_.resize(n);
}


bool Foam::seulexSparseLU::addEntry(const label d)
{
    if (n_ == 0 || col_.empty())
    {
        return false;
    }

    const label i = iperm_[d/n_];
    const label j = iperm_[d % n_];
    const label* first = col_.cdata() + rowStart_[i];
    const label* last = col_.cdata() + rowStart_[i + 1];
    const label* p = std::lower_bound(first, last, j);

    if (p == last || *p != j)
    {
        return false;
    }

    src_[p - col_.cdata()] = d;
    return true;
}


// ************************************************************************* //
