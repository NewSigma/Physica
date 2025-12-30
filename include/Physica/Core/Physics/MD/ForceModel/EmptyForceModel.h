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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "Physica/Core/Physics/MD/MDCell.h"

namespace Physica {
    template<Scalar T, unsigned int Dim> class EmptyForceModel;

    namespace Internal {
        template<class>
        struct is_empty_force_model {
            constexpr static bool value = false;
        };

        template<Scalar T, unsigned int Dim>
        struct is_empty_force_model<EmptyForceModel<T, Dim>> {
            constexpr static bool value = true;
        };
    }
    /**
     * \class EmptyForceModel serves starting point for a new force model.
     */
    template<Scalar T, unsigned int Dim>
    class EmptyForceModel {
    public:
        using MDCellType = MDCell<T, Dim>;
        using LatticeMatrix = MDCellType::LatticeMatrix;
        using ForceConstMatrix = DenseSymmMatrix<T>;
    public:
        /* Operations */
        [[nodiscard]] T potentialV(const MDCellType&) const { return 0; }

        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force(const MDCellType& cell) const { return VectorND<T>(cell.getDOF(), 0); }
        template<ExecutePolicy P>
        void forceAsync(const MDCellType& cell, Vector auto& result) const;
        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force_short(const MDCellType& cell) const { return force<P>(cell); }
        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force_long(const MDCellType& cell) const { return VectorND<T>(cell.getDOF(), 0); }

        [[nodiscard]] T forceConst(const MDCellType&, [[maybe_unused]] size_t dof1, [[maybe_unused]] size_t dof2) const { return T(0); }
        [[nodiscard]] ForceConstMatrix forceConst(const MDCellType& cell) const;

        [[nodiscard]] LatticeMatrix virial(const MDCellType&) const { return LatticeMatrix(Dim, Dim, 0); }
    };

    template<Scalar T, unsigned int Dim>
    template<ExecutePolicy P>
    void EmptyForceModel<T, Dim>::forceAsync(
            [[maybe_unused]] const MDCellType& cell, Vector auto& result) const {
        assert(result.getLength() == cell.getDOF() && "[Error]: Array length does not match");
        result = T(0);
    }

    template<Scalar T, unsigned int Dim>
    auto EmptyForceModel<T, Dim>::forceConst(const MDCellType& cell) const -> ForceConstMatrix {
        return ForceConstMatrix(cell.getDOF(), T(0));
    }
}

namespace Physica {
    template<Scalar T, unsigned int Dim>
    class Traits<EmptyForceModel<T, Dim>> {
    public:
        constexpr static bool IsPeriodBoundary = true;
        constexpr static bool IsContractable = false;
    };
}
