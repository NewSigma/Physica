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

namespace Physica::Core {
    template<Matrix T>
    class LValueFlatten<T> : public LValueVector<LValueFlatten<T>> {
        using This = LValueFlatten<T>;

        const T& mat;
    public:
        using Base = LValueVector<This>;
        using typename Base::ScalarType;
    protected:
        using PtrTy = typename ScalarType::PtrTy;
        using ConstPtrTy = typename ScalarType::ConstPtrTy;
    public:
        LValueFlatten(const LValueMatrix<T>& mat_) : mat(mat_.getDerived()) {}
        LValueFlatten(const This&) = delete;
        LValueFlatten(This&&) noexcept = delete;
        ~LValueFlatten() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operators */
        [[nodiscard]] ScalarType& operator[](size_t index) { return *data_ptr(index); }
        [[nodiscard]] const ScalarType& operator[](size_t index) const { return *data_ptr(index); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow() * mat.getCol(); }
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(size_t index);
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(size_t index) const;
    };

    template<Matrix T>
    __host__ __device__ inline typename LValueFlatten<T>::PtrTy
    LValueFlatten<T>::data_ptr(size_t index) {
        return const_cast<ScalarType*>(const_cast<const This&>(*this).data_ptr(index));
    }

    template<Matrix T>
    __host__ __device__ inline typename LValueFlatten<T>::ConstPtrTy
    LValueFlatten<T>::data_ptr(size_t index) const {
        const size_t major = index / mat.getMaxMinor();
        const size_t minor = index % mat.getMaxMinor();
        const size_t row = MatrixOption::rowFromMajorMinor<T>(major, minor);
        const size_t col = MatrixOption::colFromMajorMinor<T>(major, minor);
        return mat.data_ptr(row, col);
    }
}

namespace Physica {
    template<Core::Matrix T>
    class Traits<Core::LValueFlatten<T>> {
    public:
        using ScalarType = typename T::ScalarType;
        constexpr static size_t SizeAtCompile = T::RowAtCompile * T::ColAtCompile;

        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };
}
