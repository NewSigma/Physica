/*
 * Copyright 2024 Weibo He.
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
    template<Scalar T, class ReprType>
    HubbardMatrix<T, ReprType>::HubbardMatrix(ModelBase hubbard, ReprType repr_)
            : ModelBase(std::move(hubbard))
            , repr(std::move(repr_))
            , planProvider(NumSite, PlanFlag::Estimate) {
        assert(ModelBase::getNumSuperCellSite() == NumSite && "[Error]: Inconsistent site number");
    }

    template<Scalar T, class ReprType>
    template<class Functor>
    void HubbardMatrix<T, ReprType>::forNeighSites(Functor func, int site) const {
        if constexpr (Dim == 1)
            func(site, (site + 1) % NumSite);
        else {
            const auto& hopTargets = ModelBase::getHopIndexArray()[site];
            for (int site1 : hopTargets)
                func(site, site1);
        }
    }

    template<Scalar T, class ReprType>
    T HubbardMatrix<T, ReprType>::calc(size_t row, size_t col) const {
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
                for (int i = 0; i < int(NumSite); ++i) {
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

    template<Scalar T, class ReprType>
    T HubbardMatrix<T, ReprType>::trace() const {
        size_t numDoubleOccupy = 0;
        for (size_t i = 0; i < getNumState(); ++i)
            numDoubleOccupy += repr[i].getNumDoubleOccupy();
        return getRepelU() * T(numDoubleOccupy);
    }

    template<Scalar T, class ReprType>
    void HubbardMatrix<T, ReprType>::swap(HubbardMatrix& __restrict obj) noexcept {
        Base::swap(obj);
        ModelBase::swap(obj);
        repr.swap(obj.repr);
        planProvider.swap(obj.planProvider);
    }

    template<Scalar T, class ReprType>
    inline typename HubbardMatrix<T, ReprType>::RealType HubbardMatrix<T, ReprType>::repelElem(StateType psi) const {
        return getRepelU() * RealType(psi.getNumDoubleOccupy());
    }

    template<Scalar T, class ReprType>
    typename HubbardMatrix<T, ReprType>::RealType
    HubbardMatrix<T, ReprType>::hoppingElem(StateType rowPsi, StateType colPsi) const {
        int count = 0;
        if constexpr (Dim == 1) {
            for (int site = 0; site < int(NumSite); ++site) {
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
        }
        else {
            ModelBase::forSiteInLattice([this, &count, rowPsi, colPsi](IndexType index) {
                const auto& dims = ModelBase::getDims();
                const int site = IndexType::toIndex1D(dims, index);
                for (unsigned int dim = 0; dim < Dim; ++dim) {
                    IndexType index1 = index;
                    index1[dim] = (index1[dim] + 1) % ModelBase::getSuperSize()[dim];

                    const int site1 = IndexType::toIndex1D(dims, index1);
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
            });
        }
        return RealType(-count) * getHoppingT();
    }
}
