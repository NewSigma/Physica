/*
 * Copyright 2024 WeiBo He.
 *
 * This file is part of Physica.
 *
 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

namespace Physica::Core {
    template<class ScalarType, class ReprType>
    Hubbard1D<ScalarType, ReprType>::Hubbard1D(
            LatticeType lattice, ReprType repr_, ScalarType hoppingT_, ScalarType repelU_)
            : Base(std::move(lattice))
            , repr(std::move(repr_))
            , hoppingT(hoppingT_)
            , repelU(repelU_) {}

    template<class ScalarType, class ReprType>
    template<class VectorType>
    Vector<ScalarType> Hubbard1D<ScalarType, ReprType>::operator*(const RValueVector<VectorType>& v) const {
        const size_t length = v.getLength();
        const auto numSite = getNumSuperCellSite();
        assert(Base::getColumn() == length && "[Error]: Dimensions do not match");
        Vector<ScalarType> result(length, 0);
        for (size_t i = 0; i < length; ++i) {
            const ScalarType hoppingElem = -v.calc(i) * hoppingT;
            const auto state = repr[i];
            int numRepel = 0;
            for (unsigned int site = 0; site < numSite; ++site) {
                const auto site1 = (site + 1) % numSite;
                Base::stateAdd(result, state.hopUp(site, site1), hoppingElem);
                Base::stateAdd(result, state.hopUp(site1, site), hoppingElem);
                Base::stateAdd(result, state.hopDown(site, site1), hoppingElem);
                Base::stateAdd(result, state.hopDown(site1, site), hoppingElem);
                numRepel += state.isUpOccupy(site) && state.isDownOccupy(site);
            }
            result[i] += v.calc(i) * repelU * ScalarType(numRepel);
        }
        return result;
    }

    template<class ScalarType, class ReprType>
    ScalarType Hubbard1D<ScalarType, ReprType>::calc(size_t row, size_t col) const {
        const auto colState = repr[col];
        const auto numSite = getNumSuperCellSite();
        if (row == col) {
            int numRepel = 0;
            for (unsigned int site = 0; site < numSite; ++site)
                numRepel += colState.isUpOccupy(site) && colState.isDownOccupy(site);
            return repelU * ScalarType(numRepel);
        }

        const auto rowState = repr[row];
        for (unsigned int site = 0; site < numSite; ++site) {
            const auto site1 = (site + 1) % numSite;
            if (rowState == colState.hopUp(site, site1)
             || rowState == colState.hopUp(site1, site)
             || rowState == colState.hopDown(site, site1)
             || rowState == colState.hopDown(site1, site))
                return -hoppingT;
        }
        return ScalarType(0);
    }

    template<class ScalarType, class ReprType>
    void Hubbard1D<ScalarType, ReprType>::swap(Hubbard1D& __restrict obj) noexcept {
        Base::swap(obj);
        repr.swap(obj.repr);
        hoppingT.swap(obj.hoppingT);
        repelU.swap(obj.repelU);
    }
}
