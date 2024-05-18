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
            LatticeType lattice, ReprType repr_, RealType hoppingT_, RealType repelU_)
            : Base(std::move(lattice))
            , repr(std::move(repr_))
            , hoppingT(hoppingT_)
            , repelU(repelU_)
            , planProvider(getNumSuperCellSite(), PlanFlag::Estimate) {}

    template<class ScalarType, class ReprType>
    template<class VectorType>
    Vector<ScalarType> Hubbard1D<ScalarType, ReprType>::operator*(const RValueVector<VectorType>& v) const {
        if constexpr (IsTransInvariant)
            return static_cast<const Base&>(*this) * v;
        else {
            const size_t length = v.getLength();
            const auto numSite = getNumSuperCellSite();
            assert(Base::getColumn() == length && "[Error]: Dimensions do not match");
            Vector<ScalarType> result(length, 0);
            for (size_t i = 0; i < length; ++i) {
                const ScalarType hop = -v.calc(i) * hoppingT;
                const auto state = repr[i];
                int numRepel = 0;
                for (unsigned int site = 0; site < numSite; ++site) {
                    const auto site1 = (site + 1) % numSite;
                    Base::stateAdd(result, state.hopUp(site, site1), hop);
                    Base::stateAdd(result, state.hopUp(site1, site), hop);
                    Base::stateAdd(result, state.hopDown(site, site1), hop);
                    Base::stateAdd(result, state.hopDown(site1, site), hop);
                    numRepel += state.isUpOccupy(site) && state.isDownOccupy(site);
                }
                result[i] += v.calc(i) * (repelU * RealType(numRepel));
            }
            return result;
        }
    }

    template<class ScalarType, class ReprType>
    ScalarType Hubbard1D<ScalarType, ReprType>::calc(size_t row, size_t col) const {
        if constexpr (IsTransInvariant) {
            const int numSite = getNumSuperCellSite();
            const auto& periods = repr.getPeriods();
            const auto psi1 = repr[row];
            if (row == col) {
                const bool flag = (repr.getKIndex() == 0) || (periods[row] == numSite);
                return flag ? repelElem(psi1) : RealType(0);
            }

            auto fft = FFTType::makeEmptyFFT(numSite);
            {
                auto& rSpace = fft.getRSpace();
                auto psi2 = repr[col];
                for (int i = 0; i < numSite; ++i) {
                    RealType elem = 0;
                    if (psi1 != psi2)
                        elem = hoppingElem(psi1, psi2);
                    rSpace[i] = elem;
                    psi2 <<= 1;
                }
            }
            FFTType::transform(planProvider, fft);
            const RealType normalizer = sqrt(RealType(periods[row] * periods[col])) / RealType(numSite);
            return fft.getKSpace()[repr.getReducedK()] * normalizer;
        }
        else {
            const auto colState = repr[col];
            if (row == col)
                return repelElem(colState);
            return hoppingElem(repr[row], colState);
        }
    }

    template<class ScalarType, class ReprType>
    void Hubbard1D<ScalarType, ReprType>::swap(Hubbard1D& __restrict obj) noexcept {
        Base::swap(obj);
        repr.swap(obj.repr);
        hoppingT.swap(obj.hoppingT);
        repelU.swap(obj.repelU);
        planProvider.swap(obj.planProvider);
    }

    template<class ScalarType, class ReprType>
    typename Hubbard1D<ScalarType, ReprType>::RealType
    Hubbard1D<ScalarType, ReprType>::repelElem(StateType psi) const {
        const auto numSite = getNumSuperCellSite();
        int numRepel = 0;
        for (unsigned int site = 0; site < numSite; ++site)
            numRepel += psi.isUpOccupy(site) && psi.isDownOccupy(site);
        return repelU * RealType(numRepel);
    }

    template<class ScalarType, class ReprType>
    typename Hubbard1D<ScalarType, ReprType>::RealType
    Hubbard1D<ScalarType, ReprType>::hoppingElem(StateType rowPsi, StateType colPsi) const {
        const auto numSite = getNumSuperCellSite();
        for (unsigned int site = 0; site < numSite; ++site) {
            const auto site1 = (site + 1) % numSite;
            if (rowPsi == colPsi.hopUp(site, site1)
             || rowPsi == colPsi.hopUp(site1, site)
             || rowPsi == colPsi.hopDown(site, site1)
             || rowPsi == colPsi.hopDown(site1, site))
                return -hoppingT;
        }
        return RealType(0);
    }
}
