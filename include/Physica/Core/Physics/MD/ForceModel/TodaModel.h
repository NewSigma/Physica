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

#include "EmptyForceModel.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] T. Hatano, Phys. Rev. E 59, R1(R) (1999); https://doi.org/10.1103/PhysRevE.59.R1
     */
    template<class ScalarType, bool IsFixedBoundary, unsigned int Dim>
    class TodaModel {
        static_assert(Dim == 1, "[Error]: TodaModel must be 1-dimensional");
        using MDCellType = MDCell<ScalarType, Dim>;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using ForceConstMatrix = typename EmptyForceModel<ScalarType, Dim>::ForceConstMatrix;

        ScalarType springLength;
    public:
        TodaModel(ScalarType springLength_) : springLength(std::move(springLength_)) {}
        TodaModel(const TodaModel&) = default;
        TodaModel(TodaModel&&) noexcept = default;
        ~TodaModel() = default;
        /* Operators */
        TodaModel& operator=(TodaModel obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] ScalarType potentialEnergy(const MDCellType& cell) const;

        template<class Executor, bool IsSmallCell = false>
        [[nodiscard]] Vector<ScalarType> force(const MDCellType& cell) const;
        template<class VectorType, class Executor, bool IsSmallCell = false>
        void forceAsync(const MDCellType& cell, ContinuousVector<VectorType>& result) const;
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force_short(const MDCellType& cell) const { return force<Executor>(cell); }
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force_long(const MDCellType& cell) const { return Vector<ScalarType>(cell.getDOF(), 0); }
        
        [[nodiscard]] ScalarType forceConst(const MDCellType& cell, size_t dof1, size_t dof2) const;
        [[nodiscard]] ForceConstMatrix forceConst(const MDCellType& cell) const;

        [[nodiscard]] LatticeMatrix virial(const MDCellType& cell) const;
        void swap(TodaModel& obj) noexcept;
    };

    template<class ScalarType, bool IsFixedBoundary, unsigned int Dim>
    ScalarType TodaModel<ScalarType, IsFixedBoundary, Dim>::potentialEnergy(const MDCellType& cell) const {
        const size_t numParticle = cell.getNumParticle();
        const auto& pos = cell.getPos();
        ScalarType energy = 0;
        if constexpr (IsFixedBoundary) {
            /* First */ {
                const ScalarType delta = pos(0, 0) - springLength;
                energy = delta + exp(-delta);
            }
            for (size_t i = 0; i < numParticle - 1; ++i) {
                const ScalarType delta = pos(i + 1, 0) - pos(i, 0) - springLength;
                energy += delta + exp(-delta);
            }
            /* Last */ {
                const ScalarType delta = cell.getLattice()(0, 0) - pos(numParticle - 1, 0) - springLength;
                energy += delta + exp(-delta);
            }
        }
        else {
            for (size_t i = 0; i < numParticle; ++i) {
                const size_t i1 = (i + 1) % numParticle;
                const ScalarType delta = cell.minDistVector(i, i1).norm() - springLength;
                energy += delta + exp(-delta);
            }
        }
        return energy;
    }

    template<class ScalarType, bool IsFixedBoundary, unsigned int Dim>
    template<class Executor, bool IsSmallCell>
    Vector<ScalarType> TodaModel<ScalarType, IsFixedBoundary, Dim>::force(const MDCellType& cell) const {
        const size_t numParticle = cell.getNumParticle();
        Vector<ScalarType> result(numParticle);
        forceAsync<Vector<ScalarType>, Executor, IsSmallCell>(cell, result);
        return result;
    }

    template<class ScalarType, bool IsFixedBoundary, unsigned int Dim>
    template<class VectorType, class Executor, bool IsSmallCell>
    void TodaModel<ScalarType, IsFixedBoundary, Dim>::forceAsync(
            const MDCellType& cell, ContinuousVector<VectorType>& result) const {
        const size_t numParticle = cell.getNumParticle();
        const auto& pos = cell.getPos();
        result = ScalarType(0);
        if constexpr (IsFixedBoundary) {
            /* First */ {
                const ScalarType delta = pos(0, 0) - springLength;
                result[0] = ScalarType(-1.0) + exp(-delta);
            }
            for (size_t i = 0; i < numParticle - 1; ++i) {
                const ScalarType delta = pos(i + 1, 0) - pos(i, 0) - springLength;
                const ScalarType f = exp(-delta);
                result[i] -= f;
                result[i + 1] += f;
            }
            /* Last */ {
                const ScalarType delta = cell.getLattice()(0, 0) - pos(numParticle - 1, 0) - springLength;
                result[numParticle - 1] += ScalarType(1.0) - exp(-delta);
            }
        }
        else {
            for (size_t i = 0; i < numParticle; ++i) {
                const size_t i1 = (i + 1) % numParticle;
                const ScalarType delta = cell.minDistVector(i, i1).norm() - springLength;
                const ScalarType f = exp(-delta);
                result[i] -= f;
                result[i1] += f;
            }
        }
    }

    template<class ScalarType, bool IsFixedBoundary, unsigned int Dim>
    ScalarType TodaModel<ScalarType, IsFixedBoundary, Dim>::forceConst(const MDCellType& cell, size_t dof1, size_t dof2) const {
        const bool isNeighbor = (dof1 + 1 == dof2) || (dof2 + 1 == dof1);
        if (isNeighbor) {
            const ScalarType delta = cell.minDistVector(dof1, dof2).norm();
            return exp(-delta);
        }
        else {
            const size_t dof3 = dof1 == 0 ? cell.getNumParticle() : dof1;
            const ScalarType delta1 = cell.minDistVector(dof1, dof1 + 1).norm();
            const ScalarType delta2 = cell.minDistVector(dof1, dof3 - 1).norm();
            return exp(-delta1) + exp(-delta2);
        }
        return 0;
    }

    template<class ScalarType, bool IsFixedBoundary, unsigned int Dim>
    typename TodaModel<ScalarType, IsFixedBoundary, Dim>::ForceConstMatrix
    TodaModel<ScalarType, IsFixedBoundary, Dim>::forceConst(const MDCellType& cell) const {
        const size_t dof = cell.getDOF();
        ForceConstMatrix result(dof, ScalarType(0));
        for (size_t i = 0; i < dof; ++i) {
            const size_t i1 = (i + 1) % dof;
            const ScalarType delta = cell.minDistVector(i, i1).norm();
            const ScalarType temp = exp(-delta);
            result(i, i1) = temp;
            result(i, i) += temp;
            result(i1, i1) += temp;
        }
        return result;
    }

    template<class ScalarType, bool IsFixedBoundary, unsigned int Dim>
    typename TodaModel<ScalarType, IsFixedBoundary, Dim>::LatticeMatrix
    TodaModel<ScalarType, IsFixedBoundary, Dim>::virial(const MDCellType& cell) const {
        const size_t numParticle = cell.getNumParticle();
        const auto& pos = cell.getPos();
        ScalarType result = 0;
        if constexpr (IsFixedBoundary) {
            /* First*/ {
                const ScalarType r = pos(0, 0);
                const ScalarType delta = r - springLength;
                const ScalarType f = exp(-delta) - ScalarType(1.0);
                result += r * f;
            }
            for (size_t i = 0; i < numParticle - 1; ++i) {
                const ScalarType r = pos(i + 1, 0) - pos(i, 0);
                const ScalarType delta = r - springLength;
                const ScalarType f = exp(-delta) - ScalarType(1.0);
                result += r * f;
            }
            /* Last */ {
                const ScalarType r = cell.getLattice()(0, 0) - pos(numParticle - 1, 0);
                const ScalarType delta = r - springLength;
                const ScalarType f = exp(-delta) - ScalarType(1.0);
                result += r * f;
            }
        }
        else {
            for (size_t i = 0; i < numParticle; ++i) {
                const size_t i1 = (i + 1) % numParticle;
                const ScalarType r = cell.minDistVector(i, i1).norm();
                const ScalarType delta = r - springLength;
                const ScalarType f = exp(-delta) - ScalarType(1.0);
                result += r * f;
            }
        }
        result /= (ScalarType(Dim) * cell.getVolume());
        return LatticeMatrix{result};
    }

    template<class ScalarType, bool IsFixedBoundary, unsigned int Dim>
    void TodaModel<ScalarType, IsFixedBoundary, Dim>::swap(TodaModel& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        springLength.swap(obj.springLength);
    }
}
