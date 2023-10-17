/*
 * Copyright 2023 WeiBo He.
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

#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"

namespace Physica::Core {
    template<class ScalarType, class PosScalarType, unsigned int Dim>
    class FPUModel {
        static_assert(Dim == 1, "[Error]: FPUModel must be 1-dimensional");
        using MDCellType = MDCell<ScalarType, PosScalarType, Dim>;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        constexpr static double SpringLength = 1.0;
    public:
        template<class Executor, bool IsSmallCell = false>
        [[nodiscard]] Vector<ScalarType> force(const MDCellType& cell) const;
        template<class VectorType, class Executor, bool IsSmallCell>
        void forceAsync(const MDCellType& cell, ContinuousVector<VectorType>& result) const;
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force_short(const MDCellType& cell) const { return force<Executor>(cell); }
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force_long(const MDCellType& cell) const { return Vector<ScalarType>(cell.getDOF(), 0); }
        [[nodiscard]] ScalarType potentialEnergy(const MDCellType& cell) const;
        [[nodiscard]] LatticeMatrix virial(const MDCellType& cell) const;
    };

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class Executor, bool IsSmallCell>
    Vector<ScalarType> FPUModel<ScalarType, PosScalarType, Dim>::force(const MDCellType& cell) const {
        Vector<ScalarType> result(cell.getNumParticle(), 0);
        forceAsync<Vector<ScalarType>, Executor, IsSmallCell>(cell, result);
        return result;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class VectorType, class Executor, bool IsSmallCell>
    void FPUModel<ScalarType, PosScalarType, Dim>::forceAsync(const MDCellType& cell, ContinuousVector<VectorType>& result) const {
        static_assert(std::is_same<Executor, SequentialExecutor>::value, "[Error]: Parallelization not implemented");
        const size_t numParticle = cell.getNumParticle();
        const auto& pos = cell.getPos();
        /* First */ {
            const ScalarType delta = pos(0, 0) - SpringLength;
            const ScalarType delta2 = square(delta);
            result[0] = -delta - delta * delta2;
        }
        for (size_t i = 0; i < numParticle - 1; ++i) {
            const ScalarType delta = pos(i + 1, 0) - pos(i, 0) - SpringLength;
            const ScalarType delta2 = square(delta);
            const ScalarType f = delta + delta * delta2;
            result[i] += f;
            result[i + 1] -= f;
        }
        /* Last */ {
            const ScalarType delta = cell.getLattice()(0, 0) - pos(numParticle - 1, 0) - SpringLength;
            const ScalarType delta2 = square(delta);
            result[numParticle - 1] += delta + delta * delta2;
        }
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    ScalarType FPUModel<ScalarType, PosScalarType, Dim>::potentialEnergy(const MDCellType& cell) const {
        const size_t numParticle = cell.getNumParticle();
        const auto& pos = cell.getPos();
        ScalarType energy;
        /* First */ {
            const ScalarType delta = pos(0, 0) - SpringLength;
            const ScalarType delta2 = square(delta);
            energy = delta2 * 0.5 + square(delta2) * 0.25;
        }
        for (size_t i = 0; i < numParticle - 1; ++i) {
            const ScalarType delta = pos(i + 1, 0) - pos(i, 0) - SpringLength;
            const ScalarType delta2 = square(delta);
            energy += delta2 * 0.5 + square(delta2) * 0.25;
        }
        /* Last */ {
            const ScalarType delta = cell.getLattice()(0, 0) - pos(numParticle - 1, 0) - SpringLength;
            const ScalarType delta2 = square(delta);
            energy += delta2 * 0.5 + square(delta2) * 0.25;
        }
        return energy;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    typename FPUModel<ScalarType, PosScalarType, Dim>::LatticeMatrix
    FPUModel<ScalarType, PosScalarType, Dim>::virial(const MDCellType& cell) const {
        const size_t numParticle = cell.getNumParticle();
        const auto& pos = cell.getPos();
        ScalarType result = 0;
        /* First */ {
            const ScalarType r = pos(0, 0);
            const ScalarType delta = r - SpringLength;
            const ScalarType delta2 = square(delta);
            const ScalarType f = -delta - delta * delta2;
            result += r * f;
        }
        for (size_t i = 0; i < numParticle - 1; ++i) {
            const ScalarType r = pos(i + 1, 0) - pos(i, 0);
            const ScalarType delta = r - SpringLength;
            const ScalarType delta2 = square(delta);
            const ScalarType f = -delta - delta * delta2;
            result += r * f;
        }
        /* Last */ {
            const ScalarType r = cell.getLattice()(0, 0) - pos(numParticle - 1, 0);
            const ScalarType delta = r - SpringLength;
            const ScalarType delta2 = square(delta);
            const ScalarType f = -delta - delta * delta2;
            result += r * f;
        }
        result /= (ScalarType(Dim) * cell.getVolume());
        return LatticeMatrix{result};
    }
}
