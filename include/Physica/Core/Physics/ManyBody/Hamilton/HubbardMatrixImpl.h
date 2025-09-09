/*
 * Copyright 2024-2025 Weibo He.
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
    template<Scalar T, Representation Repr, BoundaryCond BC>
    HubbardMatrix<T, Repr, BC>::HubbardMatrix(T hoppingT_, Tr repelU_, Lattice lattice, Repr repr_)
            : Lattice(std::move(lattice))
            , hoppingT(hoppingT_)
            , repelU(repelU_)
            , repr(std::move(repr_)) {
        assert(Lattice::getNumSuperCellSite() == NumSite && "[Error]: Inconsistent site number");
        if constexpr (BC == BoundaryCond::TBC)
            phases = Lattice::template calcPhase<Tr>();
    }

    template<Scalar T, Representation Repr, BoundaryCond BC>
    auto HubbardMatrix<T, Repr, BC>::operator*(const Vector auto& v) const noexcept {
        using V = std::remove_cvref_t<decltype(v)>;
        return HubbardVecProd<T, Repr, BC, V>(*this, v);
    }

    template<Scalar T, Representation Repr, BoundaryCond BC>
    T HubbardMatrix<T, Repr, BC>::calc(size_t row, size_t col) const {
        if constexpr (Base::IsTransInvariant) {
            const auto& periods = repr.getPeriods();
            const auto psi1 = repr[row];
            if (row == col) {
                const bool flag = (repr.getKIndex() == 0) || (periods[row] == NumSite);
                return flag ? repelElem(psi1) : Tr(0);
            }

            auto fft = FFT1D(NumSite, PlanFlag::Estimate);
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
            fft.transform();
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

    template<Scalar T, Representation Repr, BoundaryCond BC>
    auto HubbardMatrix<T, Repr, BC>::calc_value(size_t row, size_t col) const -> Tv {
        return calc(row, col).value();
    }

    template<Scalar T, Representation Repr, BoundaryCond BC>
    T HubbardMatrix<T, Repr, BC>::trace() const {
        size_t numDoubleOccupy = 0;
        for (size_t i = 0; i < Base::getNumState(); ++i)
            numDoubleOccupy += repr[i].getNumDoubleOccupy();
        return getRepelU() * T(numDoubleOccupy);
    }

    template<Scalar T, Representation Repr, BoundaryCond BC>
    void HubbardMatrix<T, Repr, BC>::swap(HubbardMatrix& __restrict obj) noexcept {
        Base::swap(obj);
        Lattice::swap(obj);
        hoppingT.swap(obj.hoppingT);
        repelU.swap(obj.repelU);
        repr.swap(obj.repr);
        phases.swap(obj.phases);
    }

    template<Scalar T, Representation Repr, BoundaryCond BC>
    void HubbardMatrix<T, Repr, BC>::setPhaseArgs(ArgVector phaseArgs) noexcept requires(BC == BoundaryCond::TBC) {
        Lattice::setPhaseArgs(phaseArgs);
        phases = Lattice::template calcPhase<Tc>();
    }

    template<Scalar T, Representation Repr, BoundaryCond BC>
    auto HubbardMatrix<T, Repr, BC>::repelElem(StateType psi) const noexcept -> Tr {
        return getRepelU() * Tr(psi.getNumDoubleOccupy());
    }

    template<Scalar T, Representation Repr, BoundaryCond BC>
    auto HubbardMatrix<T, Repr, BC>::hoppingElem(StateType rowPsi, StateType colPsi) const noexcept -> T {
        if constexpr (BC == BoundaryCond::TBC) {
            T result = 0;
            for (int site = 0; site < int(NumSite); ++site) {
                Lattice::forNeighSites([this, &result, rowPsi, colPsi](int site, int site1) noexcept {
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
                Lattice::forNeighSites([&count, rowPsi, colPsi](int site, int site1) noexcept {
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

    template<Scalar T, Representation Repr, BoundaryCond BC>
    Vector2D<T> HubbardMatrix<T, Repr, BC>::calcBoundaryPhase(int site, int site1) const noexcept {
        Vector2D<T> result;
        if constexpr (BC == BoundaryCond::TBC) {
            const auto& boundary = Lattice::getSiteBoundaryMap();
            auto pair = std::make_pair(site, site1);
            if (boundary.contains(pair)) {
                int dim = boundary.find(pair)->second;
                result[0] = phases[dim];
                result[1] = phases[dim].conjugate();
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
