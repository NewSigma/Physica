/*
 * Copyright 2023-2024 WeiBo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "Physica/Core/Physics/MD/MDCell.h"

namespace Physica::Core {
    template<class ScalarType, unsigned int Dim> class EmptyForceModel;

    namespace Internal {
        template<class T>
        struct is_empty_force_model {
            constexpr static bool value = false;
        };

        template<class ScalarType, unsigned int Dim>
        struct is_empty_force_model<EmptyForceModel<ScalarType, Dim>> {
            constexpr static bool value = true;
        };

        template<class ScalarType, unsigned int Dim>
        class Traits<EmptyForceModel<ScalarType, Dim>> {
        public:
            constexpr static bool IsPeriodBoundary = true;
            constexpr static bool IsContractable = false;
        };
    }

    template<class ScalarType, unsigned int Dim>
    class EmptyForceModel {
    public:
        using MDCellType = MDCell<ScalarType, Dim>;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using ForceConstMatrix = DenseSymmMatrix<ScalarType>;
    public:
        /* Operations */
        [[nodiscard]] ScalarType potentialV([[maybe_unused]] const MDCellType& cell) const { return 0; }

        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force(const MDCellType& cell) const { return Vector<ScalarType>(cell.getDOF(), 0); }
        template<class VectorType, class Executor>
        void forceAsync([[maybe_unused]] const MDCellType& cell, ContinuousVector<VectorType>& result) const;
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force_short(const MDCellType& cell) const { return force<Executor>(cell); }
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force_long(const MDCellType& cell) const { return Vector<ScalarType>(cell.getDOF(), 0); }

        [[nodiscard]] ScalarType forceConst([[maybe_unused]] const MDCellType& cell, [[maybe_unused]] size_t dof1, [[maybe_unused]] size_t dof2) const { return ScalarType(0); }
        [[nodiscard]] ForceConstMatrix forceConst([[maybe_unused]] const MDCellType& cell) const;

        [[nodiscard]] LatticeMatrix virial([[maybe_unused]] const MDCellType& cell) const { return LatticeMatrix(Dim, Dim, 0); }
    };

    template<class ScalarType, unsigned int Dim>
    template<class VectorType, class Executor>
    void EmptyForceModel<ScalarType, Dim>::forceAsync(
            [[maybe_unused]] const MDCellType& cell, ContinuousVector<VectorType>& result) const {
        assert(result.getLength() == cell.getDOF() && "[Error]: Array length does not match");
        result = ScalarType(0);
    }

    template<class ScalarType, unsigned int Dim>
    typename EmptyForceModel<ScalarType, Dim>::ForceConstMatrix
    EmptyForceModel<ScalarType, Dim>::forceConst([[maybe_unused]] const MDCellType& cell) const {
        return ForceConstMatrix(cell.getDOF(), ScalarType(0));
    }
}
