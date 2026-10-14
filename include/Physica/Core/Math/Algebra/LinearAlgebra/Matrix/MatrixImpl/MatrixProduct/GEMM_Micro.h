/*
 * Copyright 2026 Weibo He.
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

#include "GEMM.h"

namespace Physica::Internal {
    template<Scalar T, int LhsMajor, int RhsMajor>
    PHYSICA_API void gemm2x2(const T* lhs, const T* rhs, T* target) noexcept;
    template<Scalar T, int LhsMajor, int RhsMajor>
    PHYSICA_API void gemm2x2_add(const T* lhs, const T* rhs, T* target) noexcept;

    extern template void gemm2x2<float64, MatrixMajor::Col, MatrixMajor::Col>(const float64*, const float64*, float64*) noexcept;
    extern template void gemm2x2<float64, MatrixMajor::Row, MatrixMajor::Col>(const float64*, const float64*, float64*) noexcept;
    extern template void gemm2x2<float64, MatrixMajor::Col, MatrixMajor::Row>(const float64*, const float64*, float64*) noexcept;
    extern template void gemm2x2<float64, MatrixMajor::Row, MatrixMajor::Row>(const float64*, const float64*, float64*) noexcept;

    extern template void gemm2x2<float32, MatrixMajor::Col, MatrixMajor::Col>(const float32*, const float32*, float32*) noexcept;
    extern template void gemm2x2<float32, MatrixMajor::Row, MatrixMajor::Col>(const float32*, const float32*, float32*) noexcept;
    extern template void gemm2x2<float32, MatrixMajor::Col, MatrixMajor::Row>(const float32*, const float32*, float32*) noexcept;
    extern template void gemm2x2<float32, MatrixMajor::Row, MatrixMajor::Row>(const float32*, const float32*, float32*) noexcept;

    extern template void gemm2x2_add<float64, MatrixMajor::Col, MatrixMajor::Col>(const float64*, const float64*, float64*) noexcept;
    extern template void gemm2x2_add<float64, MatrixMajor::Row, MatrixMajor::Col>(const float64*, const float64*, float64*) noexcept;
    extern template void gemm2x2_add<float64, MatrixMajor::Col, MatrixMajor::Row>(const float64*, const float64*, float64*) noexcept;
    extern template void gemm2x2_add<float64, MatrixMajor::Row, MatrixMajor::Row>(const float64*, const float64*, float64*) noexcept;

    extern template void gemm2x2_add<float32, MatrixMajor::Col, MatrixMajor::Col>(const float32*, const float32*, float32*) noexcept;
    extern template void gemm2x2_add<float32, MatrixMajor::Row, MatrixMajor::Col>(const float32*, const float32*, float32*) noexcept;
    extern template void gemm2x2_add<float32, MatrixMajor::Col, MatrixMajor::Row>(const float32*, const float32*, float32*) noexcept;
    extern template void gemm2x2_add<float32, MatrixMajor::Row, MatrixMajor::Row>(const float32*, const float32*, float32*) noexcept;
}

namespace Physica {
    template<Matrix M1, Matrix M2>
    void GEMM<M1, M2>::assign_micro(Matrix auto&& target) const noexcept {
        target.assert_assign(*this);
        using Target = std::remove_cvref_t<decltype(target)>;
        using LHS = std::remove_cvref_t<M1>;
        using RHS = std::remove_cvref_t<M2>;
        static_assert(useMicro(target));

        constexpr auto flip = [](int major) consteval static noexcept {
            return major == MatrixMajor::Row ? MatrixMajor::Col : MatrixMajor::Row;
        };
        constexpr bool TargetIsRow = Target::isRowMatrix();
        constexpr int MajorLHS = LHS::getMajor();
        constexpr int MajorRHS = RHS::getMajor();
        constexpr int KernelMajorLHS = TargetIsRow ? flip(MajorRHS) : MajorLHS;
        constexpr int KernelMajorRHS = TargetIsRow ? flip(MajorLHS) : MajorRHS;

        const auto* lhsData = getLHS().data_handle();
        const auto* rhsData = getRHS().data_handle();
        auto* targetData = target.data_handle();
        const auto* kernelLHS = TargetIsRow ? rhsData : lhsData;
        const auto* kernelRHS = TargetIsRow ? lhsData : rhsData;
        Internal::gemm2x2<T, KernelMajorLHS, KernelMajorRHS>(kernelLHS, kernelRHS, targetData);
    }

    template<Matrix M1, Matrix M2>
    void GEMM<M1, M2>::assign_add_micro(Matrix auto& target) const noexcept {
        using Target = std::remove_cvref_t<decltype(target)>;
        using LHS = std::remove_cvref_t<M1>;
        using RHS = std::remove_cvref_t<M2>;
        static_assert(useMicro(target));

        constexpr auto flip = [](int major) consteval static noexcept {
            return major == MatrixMajor::Row ? MatrixMajor::Col : MatrixMajor::Row;
        };
        constexpr bool TargetIsRow = Target::isRowMatrix();
        constexpr int MajorLHS = LHS::getMajor();
        constexpr int MajorRHS = RHS::getMajor();
        constexpr int KernelMajorLHS = TargetIsRow ? flip(MajorRHS) : MajorLHS;
        constexpr int KernelMajorRHS = TargetIsRow ? flip(MajorLHS) : MajorRHS;

        const auto* lhsData = getLHS().data_handle();
        const auto* rhsData = getRHS().data_handle();
        auto* targetData = target.data_handle();
        const auto* kernelLHS = TargetIsRow ? rhsData : lhsData;
        const auto* kernelRHS = TargetIsRow ? lhsData : rhsData;
        Internal::gemm2x2_add<T, KernelMajorLHS, KernelMajorRHS>(kernelLHS, kernelRHS, targetData);
    }

    template<Matrix M1, Matrix M2>
    consteval bool GEMM<M1, M2>::useMicro(const Matrix auto& target) noexcept {
        using Target = std::remove_cvref_t<decltype(target)>;
        using LHS = std::remove_cvref_t<M1>;
        using RHS = std::remove_cvref_t<M2>;
        if constexpr (Base::isDiffable() || Base::isComplex())
            return false;
        else {
            constexpr bool IsCompact = Target::isCompact() && LHS::isCompact() && RHS::isCompact();
            constexpr auto isFixedSize = [](int dim) consteval static noexcept {
                return LHS::getRowAtCompile() == dim
                    && LHS::getColAtCompile() == dim
                    && RHS::getRowAtCompile() == dim
                    && RHS::getColAtCompile() == dim;
            };
            return IsCompact && isFixedSize(2) && BestPacket<T, 4>::Size == 4;
        }
    }
}
