/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 Nanegate
     \\/     M anipulation  |
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::seulexSparseLU::factor
(
    const scalarSquareMatrix& J,
    const scalar r
)
{
    if (col_.empty())
    {
        return false;
    }

    const scalar* Jv = J.v();
    forAll(val_, p)
    {
        val_[p] = src_[p] >= 0 ? -Jv[src_[p]] : 0;
    }

    // Pivot scale: the diagonal before elimination (the rows mix units,
    // K against kmol/m^3, so a row maximum is no scale)
    for (label i = 0; i < n_; ++i)
    {
        val_[diag_[i]] += r;
        scale_[i] = mag(val_[diag_[i]]) + mag(r);
    }

    scalar* v = val_.data();
    const label* dst = updDst_.cdata();
    const label* src = updSrc_.cdata();

    forAll(elimP_, e)
    {
        const label p = elimP_[e];
        const label cnt = elimCnt_[e];
        if (v[p] != 0)
        {
            const scalar piv = v[elimPiv_[e]];
            if (!(mag(piv) > 1e-10*scale_[col_[p]]))
            {
                return false;
            }
            const scalar l = v[p]/piv;
            v[p] = l;
            for (label c = 0; c < cnt; ++c)
            {
                v[dst[c]] -= l*v[src[c]];
            }
        }
        dst += cnt;
        src += cnt;
    }

    for (label i = 0; i < n_; ++i)
    {
        if (!(mag(v[diag_[i]]) > 1e-10*scale_[i]))
        {
            return false;
        }
    }

    return true;
}


inline void Foam::seulexSparseLU::solve(scalarField& b)
{
    scalar* w = work_.data();

    // Forward substitution from the first non-zero entry of the permuted
    // right-hand side: a leading zero block stays zero
    label firstNz = n_;
    for (label i = 0; i < n_; ++i)
    {
        w[i] = b[perm_[i]];
        if (firstNz == n_ && w[i] != 0)
        {
            firstNz = i;
        }
    }
    for (label i = firstNz; i < n_; ++i)
    {
        scalar s = w[i];
        for (label p = rowStart_[i]; p < diag_[i]; ++p)
        {
            s -= val_[p]*w[col_[p]];
        }
        w[i] = s;
    }

    // Backward substitution
    for (label i = n_ - 1; i >= 0; --i)
    {
        scalar s = w[i];
        for (label p = upper_[i]; p < rowStart_[i + 1]; ++p)
        {
            s -= val_[p]*w[col_[p]];
        }
        w[i] = s/val_[diag_[i]];
    }

    for (label i = 0; i < n_; ++i)
    {
        b[perm_[i]] = w[i];
    }
}


// ************************************************************************* //
