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

#include "HubbardMatrix.h"

namespace Physica {
    template<Scalar T, Representation U, BoundaryCond BC>
    HubbardMatrix<T, U, BC>::HubbardMatrix(ModelBase hubbard, U repr_)
            : ModelBase(std::move(hubbard))
            , repr(std::move(repr_))
            , planProvider(NumSite, PlanFlag::Estimate) {
        assert(ModelBase::getNumSuperCellSite() == NumSite && "[Error]: Inconsistent site number");
    }

    template<Scalar T, Representation U, BoundaryCond BC>
    HubbardMatrix<T, U, BC>::HubbardMatrix(ModelBase hubbard, U repr_, DenseVector<Tr, Dim> phaseArgs_) requires(BC == BoundaryCond::TBC)
            : HubbardMatrix(std::move(hubbard), std::move(repr_)) {
        phaseArgs = std::move(phaseArgs_);
    }

    template<Scalar T, Representation U, BoundaryCond BC>
    auto HubbardMatrix<T, U, BC>::operator*(const Vector auto& v) const noexcept {
        using V = std::remove_cvref_t<decltype(v)>;
        return HubbardVecProd<T, U, BC, V>(*this, v);
    }

    template<Scalar T, Representation U, BoundaryCond BC>
    T HubbardMatrix<T, U, BC>::calc(size_t row, size_t col) const {
        if constexpr (Base::IsTransInvariant) {
            const auto& periods = repr.getPeriods();
            const auto psi1 = repr[row];
            if (row == col) {
                const bool flag = (repr.getKIndex() == 0) || (periods[row] == NumSite);
                return flag ? repelElem(psi1) : Tr(0);
            }

            auto fft = FFT1D::makeEmptyFFT(NumSite);
            {
                auto& rSpace = fft.getRSpace();
                auto psi2 = repr[col];
                int sign = 1;
                for (int i = 0; i < int(NumSite); ++i) {
                    Tr elem = 0;
                    if (psi1 != psi2)
                        elem = hoppingElem(psi1, psi2).real() * Tr(sign);
                    rSpace[i] = elem;

                    sign *= psi2.lShiftSign();
                    psi2 <<= 1;
                }
            }
            FFT1D::transform(planProvider, fft);
            const Tr normalizer = sqrt(Tr(periods[row] * periods[col])) / Tr(NumSite);
            return fft.getKSpace()[repr.getReducedK()] * normalizer;
        }
        else {
            const auto colState = repr[col];
            if (row == col)
                return repelElem(colState);
            return hoppingElem(repr[row], colState);
        }
    }

    template<Scalar T, Representation U, BoundaryCond BC>
    T HubbardMatrix<T, U, BC>::trace() const {
        size_t numDoubleOccupy = 0;
        for (size_t i = 0; i < Base::getNumState(); ++i)
            numDoubleOccupy += repr[i].getNumDoubleOccupy();
        return getRepelU() * T(numDoubleOccupy);
    }

    template<Scalar T, Representation U, BoundaryCond BC>
    void HubbardMatrix<T, U, BC>::swap(HubbardMatrix& __restrict obj) noexcept {
        Base::swap(obj);
        ModelBase::swap(obj);
        repr.swap(obj.repr);
        planProvider.swap(obj.planProvider);
    }

    template<Scalar T, Representation U, BoundaryCond BC>
    auto HubbardMatrix<T, U, BC>::repelElem(StateType psi) const noexcept -> Tr {
        return getRepelU() * Tr(psi.getNumDoubleOccupy());
    }

    template<Scalar T, Representation U, BoundaryCond BC>
    auto HubbardMatrix<T, U, BC>::hoppingElem(StateType rowPsi, StateType colPsi) const noexcept -> T {
        if constexpr (BC == BoundaryCond::TBC) {
            T result = 0;
            for (int site = 0; site < int(NumSite); ++site) {
                ModelBase::forNeighSites([this, &result, rowPsi, colPsi](int site, int site1) noexcept {
                    const int signUp = colPsi.hopUpSign(site, site1);
                    const int signDown = colPsi.hopDownSign(site, site1);
                    Vector2D<T> phases = calcBoundaryPhase(site, site1);

                    T n1 = phases[0] * Tr(rowPsi == colPsi.hopUp(site, site1))
                         - phases[1] * Tr(rowPsi == colPsi.hopUp(site1, site));
                    T n2 = phases[0] * Tr(rowPsi == colPsi.hopDown(site, site1))
                         - phases[1] * Tr(rowPsi == colPsi.hopDown(site1, site));
                    result += n1 * Tr(signUp) + n2 * Tr(signDown);
                }, site);
            }
            return -result * getHoppingT();
        }
        else {
            int count = 0;
            for (int site = 0; site < int(NumSite); ++site) {
                ModelBase::forNeighSites([&count, rowPsi, colPsi](int site, int site1) noexcept {
                    const int signUp = colPsi.hopUpSign(site, site1);
                    const int signDown = colPsi.hopDownSign(site, site1);
                    int n1 = (rowPsi == colPsi.hopUp(site, site1))
                           - (rowPsi == colPsi.hopUp(site1, site));
                    int n2 = (rowPsi == colPsi.hopDown(site, site1))
                           - (rowPsi == colPsi.hopDown(site1, site));
                    count += n1 * signUp + n2 * signDown;
                }, site);
            }
            return Tr(-count) * getHoppingT();
        }
    }

    template<Scalar T, Representation U, BoundaryCond BC>
    Vector2D<T> HubbardMatrix<T, U, BC>::calcBoundaryPhase(int site, int site1) const noexcept {
        Vector2D<T> result;
        if constexpr (BC == BoundaryCond::TBC) {
            const auto& map = ModelBase::getSiteBoundaryMap();
            if (map.contains(std::make_pair(site, site1))) {
                int dim = map.find(std::make_pair(site, site1))->second;
                result[0] = T::fromPhase(phaseArgs[dim]);
                result[1] = T::fromPhase(-phaseArgs[dim]);
            }
            else if (map.contains(std::make_pair(site1, site))) {
                int dim = map.find(std::make_pair(site1, site))->second;
                result[0] = T::fromPhase(-phaseArgs[dim]);
                result[1] = T::fromPhase(phaseArgs[dim]);
            }
            else
                result = Tr(1);
        }
        else {
            static_assert(BC == BoundaryCond::PBC, "[Error]: Unsupported boundary condition");
            result = Tr(1);
        }
        return result;
    }
}
