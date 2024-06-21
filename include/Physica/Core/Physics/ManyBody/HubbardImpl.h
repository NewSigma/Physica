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
    Hubbard<ScalarType, ReprType>::Hubbard(
            DimArray superSize, unsigned int numUnitCellSite, ReprType repr_, RealType hoppingT_, RealType repelU_)
            : Base(std::move(superSize), numUnitCellSite)
            , repr(std::move(repr_))
            , hoppingT(hoppingT_)
            , repelU(repelU_)
            , planProvider(NumSite, PlanFlag::Estimate) {}

    template<class ScalarType, class ReprType>
    template<class AnyVector>
    Vector<ScalarType> Hubbard<ScalarType, ReprType>::operator*(const RValueVector<AnyVector>& v) const {
        const size_t length = v.getLength();
        assert(Base::getColumn() == length && "[Error]: Dimensions do not match");
        Vector<ScalarType> result(length, 0);
        for (size_t i = 0; i < length; ++i) {
            if constexpr (Dim == 1)
                dotImpl1D(result, v.calc(i), i);
            else
                dotImplND(result, v.calc(i), i);
        }
        return result;
    }

    template<class ScalarType, class ReprType>
    ScalarType Hubbard<ScalarType, ReprType>::calc(size_t row, size_t col) const {
        if constexpr (IsTransInvariant) {
            const auto& periods = repr.getPeriods();
            const auto psi1 = repr[row];
            if (row == col) {
                const bool flag = (repr.getKIndex() == 0) || (periods[row] == NumSite);
                return flag ? repelElem(psi1) : RealType(0);
            }

            auto fft = FFTType::makeEmptyFFT(NumSite);
            {
                auto& rSpace = fft.getRSpace();
                auto psi2 = repr[col];
                int sign = 1;
                for (int i = 0; i < static_cast<int>(NumSite); ++i) {
                    RealType elem = 0;
                    if (psi1 != psi2)
                        elem = hoppingElem(psi1, psi2) * RealType(sign);
                    rSpace[i] = elem;

                    sign *= psi2.lShiftSign();
                    psi2 <<= 1;
                }
            }
            FFTType::transform(planProvider, fft);
            const RealType normalizer = sqrt(RealType(periods[row] * periods[col])) / RealType(NumSite);
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
    void Hubbard<ScalarType, ReprType>::swap(Hubbard& __restrict obj) noexcept {
        Base::swap(obj);
        repr.swap(obj.repr);
        hoppingT.swap(obj.hoppingT);
        repelU.swap(obj.repelU);
        planProvider.swap(obj.planProvider);
    }

    template<class ScalarType, class ReprType>
    typename Hubbard<ScalarType, ReprType>::RealType
    Hubbard<ScalarType, ReprType>::repelElem(StateType psi) const {
        int numRepel = 0;
        for (unsigned int site = 0; site < NumSite; ++site)
            numRepel += psi.isUpOccupy(site) && psi.isDownOccupy(site);
        return repelU * RealType(numRepel);
    }

    template<class ScalarType, class ReprType>
    typename Hubbard<ScalarType, ReprType>::RealType
    Hubbard<ScalarType, ReprType>::hoppingElem(StateType rowPsi, StateType colPsi) const {
        int count = 0;
        for (unsigned int site = 0; site < NumSite; ++site) {
            const auto site1 = (site + 1) % NumSite;
            const int signUp = colPsi.hopUpSign(site, site1);
            const int signDown = colPsi.hopDownSign(site, site1);
            int n1 = 0;
            n1 += rowPsi == colPsi.hopUp(site, site1);
            n1 -= rowPsi == colPsi.hopUp(site1, site);
            int n2 = 0;
            n2 += rowPsi == colPsi.hopDown(site, site1);
            n2 -= rowPsi == colPsi.hopDown(site1, site);
            count += n1 * signUp + n2 * signDown;
        }
        return RealType(-count) * hoppingT;
    }

    template<class ScalarType, class ReprType>
    void Hubbard<ScalarType, ReprType>::sumHopping(Vector<ScalarType>& target, ScalarType value, StateType psi) const noexcept {
        if (psi.isVacuum())
            return;
        const auto index = getRepr()[psi];
        target[index] += value;
    }

    template<class ScalarType, class ReprType>
    void Hubbard<ScalarType, ReprType>::sumHopping(VectorType& target, FFTType& fft, ScalarType factor, StateType psi) const {
        if (psi.isVacuum())
            return;
        const auto reducedPsi = psi.transReduce();
        auto& rSpace = fft.getRSpace();
        int sign = 1;
        for (size_t i = 0; i < NumSite; ++i) {
            rSpace[i] = RealType(reducedPsi == psi ? sign : 0);
            sign *= psi.lShiftSign();
            psi <<= 1;
        }
        FFTType::transform(planProvider, fft);
        const size_t index = repr[reducedPsi];
        target[index] += fft.getKSpace()[repr.getReducedK()] * sqrt(RealType(repr.getPeriods()[index])) * factor;
    }

    template<class ScalarType, class ReprType>
    void Hubbard<ScalarType, ReprType>::dotImpl1D(Vector<ScalarType>& result, ScalarType factor, size_t index) const {
        ScalarType hop = -factor * hoppingT;
        if constexpr (IsTransInvariant) {
            const RealType normalizer = sqrt(RealType(repr.getPeriods()[index])) / RealType(NumSite);
            hop *= normalizer;
        }

        const auto state = repr[index];
        int numRepel = 0;
        if constexpr (IsTransInvariant) {
            auto fft = FFTType::makeEmptyFFT(NumSite);
            for (unsigned int site = 0; site < NumSite; ++site) {
                const auto site1 = (site + 1) % NumSite;
                const ScalarType hopUp = hop * RealType(state.hopUpSign(site, site1));
                const ScalarType hopDown = hop * RealType(state.hopDownSign(site, site1));
                sumHopping(result, fft, hopUp, state.hopUp(site, site1));
                sumHopping(result, fft, -hopUp, state.hopUp(site1, site));
                sumHopping(result, fft, hopDown, state.hopDown(site, site1));
                sumHopping(result, fft, -hopDown, state.hopDown(site1, site));
                numRepel += state.isUpOccupy(site) && state.isDownOccupy(site);
            }
        }
        else {
            for (unsigned int site = 0; site < NumSite; ++site) {
                const auto site1 = (site + 1) % NumSite;
                const ScalarType hopUp = hop * RealType(state.hopUpSign(site, site1));
                const ScalarType hopDown = hop * RealType(state.hopDownSign(site, site1));
                sumHopping(result, hopUp, state.hopUp(site, site1));
                sumHopping(result, -hopUp, state.hopUp(site1, site));
                sumHopping(result, hopDown, state.hopDown(site, site1));
                sumHopping(result, -hopDown, state.hopDown(site1, site));
                numRepel += state.isUpOccupy(site) && state.isDownOccupy(site);
            }
        }
        result[index] += factor * (repelU * RealType(numRepel));
    }

    template<class ScalarType, class ReprType>
    void Hubbard<ScalarType, ReprType>::dotImplND(Vector<ScalarType>& result, ScalarType factor, size_t index) const {
        static_assert(!IsTransInvariant && "[Error]: Not implemented");
        const ScalarType hop = -factor * hoppingT;
        const auto state = repr[index];
        unsigned char numRepel = 0;
        Base::forSiteInLattice([this, &result, &state, &numRepel, hop](IndexType site) {
            for (unsigned int dim = 0; dim < Dim; ++dim) {
                SiteIndex site1 = site;
                site1[dim] = (site1[dim] + 1) % Base::getSuperSize()[dim];
                const ScalarType hopUp = hop * RealType(state.hopUpSign(site, site1));
                const ScalarType hopDown = hop * RealType(state.hopDownSign(site, site1));
                sumHopping(result, hopUp, state.hopUp(site, site1));
                sumHopping(result, -hopUp, state.hopUp(site1, site));
                sumHopping(result, hopDown, state.hopDown(site, site1));
                sumHopping(result, -hopDown, state.hopDown(site1, site));
            }
            numRepel += state.isUpOccupy(site) && state.isDownOccupy(site);
        });
        result[index] += factor * (repelU * RealType(numRepel));
    }
}
