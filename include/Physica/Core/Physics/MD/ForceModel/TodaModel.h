/*
 * Copyright 2023-2024 Weibo He.
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

#include "EmptyForceModel.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] T. Hatano, Phys. Rev. E 59, R1(R) (1999); https://doi.org/10.1103/PhysRevE.59.R1
     */
    template<Scalar T, bool IsPeriodBoundary>
    class TodaModel {
        using MDCellType = MDCell<T, 1>;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using ForceConstMatrix = typename EmptyForceModel<T, 1>::ForceConstMatrix;

        T springLength;
    public:
        TodaModel(T springLength_) : springLength(std::move(springLength_)) {}
        TodaModel(const TodaModel&) = default;
        TodaModel(TodaModel&&) noexcept = default;
        ~TodaModel() = default;
        /* Operators */
        TodaModel& operator=(TodaModel obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] T potentialV(const MDCellType& cell) const;

        template<class Executor>
        [[nodiscard]] VectorND<T> force(const MDCellType& cell) const;
        template<Vector V, class Executor>
        void forceAsync(const MDCellType& cell, ContinuousVector<V>& result) const;
        template<class Executor>
        [[nodiscard]] VectorND<T> force_short(const MDCellType& cell) const { return force<Executor>(cell); }
        template<class Executor>
        [[nodiscard]] VectorND<T> force_long(const MDCellType& cell) const { return VectorND<T>(cell.getDOF(), 0); }
        
        [[nodiscard]] T forceConst(const MDCellType& cell, size_t dof1, size_t dof2) const;
        [[nodiscard]] ForceConstMatrix forceConst(const MDCellType& cell) const;

        [[nodiscard]] LatticeMatrix virial(const MDCellType& cell) const;
        void swap(TodaModel& __restrict obj) noexcept;
    };

    template<Scalar T, bool IsPeriodBoundary>
    T TodaModel<T, IsPeriodBoundary>::potentialV(const MDCellType& cell) const {
        const size_t numParticle = cell.getNumParticle();
        const auto& pos = cell.getPos();
        T energy = 0;
        if constexpr (IsPeriodBoundary) {
            for (size_t i = 0; i < numParticle; ++i) {
                const size_t i1 = (i + 1) % numParticle;
                const T delta = cell.minDistVector(i, i1).norm() - springLength;
                energy += delta + exp(-delta);
            }
        }
        else {
            /* First */ {
                const T delta = pos(0, 0) - springLength;
                energy = delta + exp(-delta);
            }
            for (size_t i = 0; i < numParticle - 1; ++i) {
                const T delta = pos(i + 1, 0) - pos(i, 0) - springLength;
                energy += delta + exp(-delta);
            }
            /* Last */ {
                const T delta = cell.getLattice()(0, 0) - pos(numParticle - 1, 0) - springLength;
                energy += delta + exp(-delta);
            }
        }
        return energy;
    }

    template<Scalar T, bool IsPeriodBoundary>
    template<class Executor>
    VectorND<T> TodaModel<T, IsPeriodBoundary>::force(const MDCellType& cell) const {
        const size_t numParticle = cell.getNumParticle();
        VectorND<T> result(numParticle);
        forceAsync<VectorND<T>, Executor>(cell, result);
        return result;
    }

    template<Scalar T, bool IsPeriodBoundary>
    template<Vector V, class Executor>
    void TodaModel<T, IsPeriodBoundary>::forceAsync(
            const MDCellType& cell, ContinuousVector<V>& result) const {
        const size_t numParticle = cell.getNumParticle();
        const auto& pos = cell.getPos();
        result = T(0);
        if constexpr (IsPeriodBoundary) {
            for (size_t i = 0; i < numParticle; ++i) {
                const size_t i1 = (i + 1) % numParticle;
                const T delta = cell.minDistVector(i, i1).norm() - springLength;
                const T f = exp(-delta);
                result[i] -= f;
                result[i1] += f;
            }
        }
        else {
            /* First */ {
                const T delta = pos(0, 0) - springLength;
                result[0] = T(-1.0) + exp(-delta);
            }
            for (size_t i = 0; i < numParticle - 1; ++i) {
                const T delta = pos(i + 1, 0) - pos(i, 0) - springLength;
                const T f = exp(-delta);
                result[i] -= f;
                result[i + 1] += f;
            }
            /* Last */ {
                const T delta = cell.getLattice()(0, 0) - pos(numParticle - 1, 0) - springLength;
                result[numParticle - 1] += T(1.0) - exp(-delta);
            }
        }
    }

    template<Scalar T, bool IsPeriodBoundary>
    T TodaModel<T, IsPeriodBoundary>::forceConst(const MDCellType& cell, size_t dof1, size_t dof2) const {
        if constexpr (IsPeriodBoundary) {
            const size_t dof = cell.getDOF();
            const bool isNeighbor = ((dof1 + 1) % dof == dof2) || ((dof2 + 1) % dof == dof1);
            if (isNeighbor) {
                const T delta = cell.minDistVector(dof1, dof2).norm() - springLength;
                return -exp(-delta);
            }

            if (dof1 == dof2) {
                const size_t dof3 = (dof1 + 1) % dof;
                const size_t dof4 = (dof1 == 0 ? cell.getNumParticle() : dof1) - 1;
                const T delta1 = cell.minDistVector(dof1, dof3).norm() - springLength;
                const T delta2 = cell.minDistVector(dof1, dof4).norm() - springLength;
                return exp(-delta1) + exp(-delta2);
            }
        }
        else {
            const auto& pos = cell.getPos();
            const bool isNeighbor = (dof1 + 1 == dof2) || (dof2 + 1 == dof1);
            if (isNeighbor) {
                const T delta = abs(pos(dof1, 0) - pos(dof2, 0)) - springLength;
                return -exp(-delta);
            }

            if (dof1 == dof2) {
                const size_t numParticle = cell.getNumParticle();
                T delta1, delta2;
                if (dof1 == 0) {
                    const T temp = pos(0, 0);
                    delta1 = temp;
                    delta2 = pos(1, 0) - temp;
                }
                else if (dof1 == numParticle - 1) {
                    const T temp = pos(numParticle - 1, 0);
                    delta1 = temp - pos(numParticle - 2, 0);
                    delta2 = cell.getLattice()(0, 0) - temp;
                }
                else {
                    const T temp = pos(dof1, 0);
                    delta1 = temp - pos(dof1 - 1, 0);
                    delta2 = pos(dof1 + 1, 0) - temp;
                }
                delta1 -= springLength;
                delta2 -= springLength;
                return exp(-delta1) + exp(-delta2);
            }
        }
        return 0;
    }

    template<Scalar T, bool IsPeriodBoundary>
    typename TodaModel<T, IsPeriodBoundary>::ForceConstMatrix
    TodaModel<T, IsPeriodBoundary>::forceConst(const MDCellType& cell) const {
        const size_t dof = cell.getDOF();
        ForceConstMatrix result(dof, T(0));
        if constexpr (IsPeriodBoundary) {
            for (size_t i = 0; i < dof; ++i) {
                const size_t i1 = (i + 1) % dof;
                const T delta = cell.minDistVector(i, i1).norm() - springLength;
                const T temp = exp(-delta);
                result(i, i1) = -temp;
                result(i, i) += temp;
                result(i1, i1) += temp;
            }
        }
        else {
            const auto& pos = cell.getPos();
            const size_t dof1 = dof - 1;
            for (size_t i = 0; i < dof1; ++i) {
                const size_t i1 = i + 1;
                const T delta = pos(i1, 0) - pos(i, 0) - springLength;
                const T temp = exp(-delta);
                result(i, i1) = -temp;
                result(i, i) += temp;
                result(i1, i1) += temp;
            }
            result(0, 0) += exp(-pos(0, 0) + springLength);
            result(dof1, dof1) += exp(-(cell.getLattice()(0, 0) - pos(dof1, 0) - springLength));
        }
        return result;
    }

    template<Scalar T, bool IsPeriodBoundary>
    typename TodaModel<T, IsPeriodBoundary>::LatticeMatrix
    TodaModel<T, IsPeriodBoundary>::virial(const MDCellType& cell) const {
        const size_t numParticle = cell.getNumParticle();
        const auto& pos = cell.getPos();
        T result = 0;
        if constexpr (IsPeriodBoundary) {
            for (size_t i = 0; i < numParticle; ++i) {
                const size_t i1 = (i + 1) % numParticle;
                const T r = cell.minDistVector(i, i1).norm();
                const T delta = r - springLength;
                const T f = exp(-delta) - T(1.0);
                result += r * f;
            }
        }
        else {
            /* First*/ {
                const T r = pos(0, 0);
                const T delta = r - springLength;
                const T f = exp(-delta) - T(1.0);
                result += r * f;
            }
            for (size_t i = 0; i < numParticle - 1; ++i) {
                const T r = pos(i + 1, 0) - pos(i, 0);
                const T delta = r - springLength;
                const T f = exp(-delta) - T(1.0);
                result += r * f;
            }
            /* Last */ {
                const T r = cell.getLattice()(0, 0) - pos(numParticle - 1, 0);
                const T delta = r - springLength;
                const T f = exp(-delta) - T(1.0);
                result += r * f;
            }
        }
        result /= cell.getVolume();
        return LatticeMatrix{result};
    }

    template<Scalar T, bool IsPeriodBoundary>
    void TodaModel<T, IsPeriodBoundary>::swap(TodaModel& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        springLength.swap(obj.springLength);
    }
}

namespace Physica {
    template<class T, bool B>
    class Traits<Core::TodaModel<T, B>> {
    public:
        using ScalarType = T;
        constexpr static bool IsPeriodBoundary = B;
        constexpr static bool IsContractable = false;
    };
}
