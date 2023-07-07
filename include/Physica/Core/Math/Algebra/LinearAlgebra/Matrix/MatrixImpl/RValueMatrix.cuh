/*
 * Copyright 2022-2023 WeiBo He.
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
    template<class Derived>
    class device_obj<RValueMatrix<Derived>> : public Utils::CRTPBase<device_obj<Derived>> {
        using Base = Utils::CRTPBase<Derived>;
        using host_obj = RValueMatrix<Derived>;
    public:
        using ScalarType = typename host_obj::ScalarType;
        constexpr static int Option = host_obj::Option;
        constexpr static size_t RowAtCompile = host_obj::RowAtCompile;
        constexpr static size_t ColumnAtCompile = host_obj::ColumnAtCompile;
        constexpr static size_t MaxRowAtCompile = host_obj::MaxRowAtCompile;
        constexpr static size_t MaxColumnAtCompile = host_obj::MaxColumnAtCompile;
        constexpr static size_t SizeAtCompile = host_obj::SizeAtCompile;
        constexpr static size_t MaxSizeAtCompile = host_obj::MaxSizeAtCompile;

        constexpr static bool isColumnMatrix = host_obj::isColumnMatrix;
        constexpr static bool isRowMatrix = host_obj::isRowMatrix;
    public:
        /* Operations */
        template<class OtherDerived>
        void assignTo(device_obj<LValueMatrix<OtherDerived>>& target) const;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const { return Base::getDerived().calc(row, col); }
        [[nodiscard]] __device__ ScalarType calcFromMajorMinor(size_t row, size_t col) const;
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return Base::getDerived().getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const noexcept { return Base::getDerived().getColumn(); }
        [[nodiscard]] __host__ __device__ inline size_t getMaxMajor() const noexcept;
        [[nodiscard]] __host__ __device__ inline size_t getMaxMinor() const noexcept;
    };
}

#ifdef __CUDA_ARCH__
    #include "RValueMatrixImpl.cuh"
#endif
