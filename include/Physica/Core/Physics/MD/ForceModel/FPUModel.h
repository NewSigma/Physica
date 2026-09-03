/*
 * Copyright 2023-2026 Weibo He.
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

#include "Physica/Core/Parallel/Parallel.h"

namespace Physica {
    template<Scalar T, bool IsPeriodBoundary, unsigned int Dim>
    class FPUModel {
        static_assert(Dim == 1, "[Error]: FPUModel must be 1-dimensional");
        using This = FPUModel<T, IsPeriodBoundary, Dim>;
        using MDCellType = MDCell<T, Dim>;
        using LatticeMatrix = MDCellType::LatticeMatrix;

        T springLength;
    public:
        FPUModel(T springLength_) : springLength(std::move(springLength_)) {}
        FPUModel(const This&) = default;
        FPUModel(This&&) noexcept = default;
        ~FPUModel() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force(const MDCellType& cell) const;
        template<ExecutePolicy P>
        void forceAsync(const MDCellType& cell, Vector auto& result) const;
        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force_short(const MDCellType& cell) const { return force<P>(cell); }
        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force_long(const MDCellType& cell) const { return VectorND<T>(cell.getDOF(), 0); }
        [[nodiscard]] T potentialV(const MDCellType& cell) const;
        [[nodiscard]] LatticeMatrix virial(const MDCellType& cell) const;
        void swap(This& __restrict obj) noexcept;
    };

    template<Scalar T, bool IsPeriodBoundary, unsigned int Dim>
    template<ExecutePolicy P>
    VectorND<T> FPUModel<T, IsPeriodBoundary, Dim>::force(const MDCellType& cell) const {
        VectorND<T> result(cell.getNumParticle());
        forceAsync<P>(cell, result);
        return result;
    }

    template<Scalar T, bool IsPeriodBoundary, unsigned int Dim>
    template<ExecutePolicy P>
    void FPUModel<T, IsPeriodBoundary, Dim>::forceAsync(const MDCellType& cell, Vector auto& result) const {
        static_assert(P == Sequential, "[Error]: Parallelization not implemented");
        const size_t numParticle = cell.getNumParticle();
        const auto& pos = cell.getPos();
        result.zeros();
        if constexpr (IsPeriodBoundary) {
            for (size_t i = 0; i < numParticle; ++i) {
                const size_t i1 = (i + 1) % numParticle;
                const T delta = cell.minDistVector(i, i1).norm() - springLength;
                const T delta2 = square(delta);
                const T f = delta + delta * delta2;
                result[i] += f;
                result[i1] -= f;
            }
        }
        else {
            /* First */ {
                const T delta = pos[0, 0] - springLength;
                const T delta2 = square(delta);
                result[0] = -delta - delta * delta2;
            }
            for (size_t i = 0; i < numParticle - 1; ++i) {
                const T delta = pos(i + 1, 0) - pos(i, 0) - springLength;
                const T delta2 = square(delta);
                const T f = delta + delta * delta2;
                result[i] += f;
                result[i + 1] -= f;
            }
            /* Last */ {
                const T delta = cell.getLattice()[0, 0] - pos(numParticle - 1, 0) - springLength;
                const T delta2 = square(delta);
                result[numParticle - 1] += delta + delta * delta2;
            }
        }
    }

    template<Scalar T, bool IsPeriodBoundary, unsigned int Dim>
    T FPUModel<T, IsPeriodBoundary, Dim>::potentialV(const MDCellType& cell) const {
        const size_t numParticle = cell.getNumParticle();
        const auto& pos = cell.getPos();
        T energy = 0;
        if constexpr (IsPeriodBoundary) {
            for (size_t i = 0; i < numParticle; ++i) {
                const size_t i1 = (i + 1) % numParticle;
                const T delta = cell.minDistVector(i, i1).norm() - springLength;
                const T delta2 = square(delta);
                energy += delta2 * 0.5 + square(delta2) * 0.25;
            }
        }
        else {
            /* First */ {
                const T delta = pos[0, 0] - springLength;
                const T delta2 = square(delta);
                energy = delta2 * 0.5 + square(delta2) * 0.25;
            }
            for (size_t i = 0; i < numParticle - 1; ++i) {
                const T delta = pos(i + 1, 0) - pos(i, 0) - springLength;
                const T delta2 = square(delta);
                energy += delta2 * 0.5 + square(delta2) * 0.25;
            }
            /* Last */ {
                const T delta = cell.getLattice()[0, 0] - pos(numParticle - 1, 0) - springLength;
                const T delta2 = square(delta);
                energy += delta2 * 0.5 + square(delta2) * 0.25;
            }
        }
        return energy;
    }

    template<Scalar T, bool IsPeriodBoundary, unsigned int Dim>
    auto FPUModel<T, IsPeriodBoundary, Dim>::virial(const MDCellType& cell) const -> LatticeMatrix {
        const size_t numParticle = cell.getNumParticle();
        const auto& pos = cell.getPos();
        T result = 0;
        if constexpr (IsPeriodBoundary) {
            for (size_t i = 0; i < numParticle; ++i) {
                const size_t i1 = (i + 1) % numParticle;
                const T r = cell.minDistVector(i, i1).norm();
                const T delta = r - springLength;
                const T delta2 = square(delta);
                const T f = -delta - delta * delta2;
                result += r * f;
            }
        }
        else {
            /* First */ {
                const T r = pos[0, 0];
                const T delta = r - springLength;
                const T delta2 = square(delta);
                const T f = -delta - delta * delta2;
                result += r * f;
            }
            for (size_t i = 0; i < numParticle - 1; ++i) {
                const T r = pos(i + 1, 0) - pos(i, 0);
                const T delta = r - springLength;
                const T delta2 = square(delta);
                const T f = -delta - delta * delta2;
                result += r * f;
            }
            /* Last */ {
                const T r = cell.getLattice()[0, 0] - pos(numParticle - 1, 0);
                const T delta = r - springLength;
                const T delta2 = square(delta);
                const T f = -delta - delta * delta2;
                result += r * f;
            }
        }
        result /= cell.getVolume();
        return LatticeMatrix{result};
    }

    template<Scalar T, bool IsPeriodBoundary, unsigned int Dim>
    void FPUModel<T, IsPeriodBoundary, Dim>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        springLength.swap(obj.springLength);
    }
}
