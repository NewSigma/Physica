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
    class device_obj<ContinuousFlatten<T>> : public device_obj<ContinuousVector<ContinuousFlatten<T>>> {
        using host_obj = ContinuousFlatten<T>;
        using This = device_obj<host_obj>;

        device_obj<T>& mat;
    public:
        using Base = device_obj<ContinuousVector<host_obj>>;
        using typename Base::ScalarType;
    protected:
        using PtrTy = typename ScalarType::PtrTy;
        using ConstPtrTy = typename ScalarType::ConstPtrTy;
    public:
        device_obj(device_obj<ContinuousMatrix<T>>& mat_) : mat(mat_.getDerived()) {}
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operators */
        [[nodiscard]] __device__ ScalarType& operator[](size_t index) { return *data_ptr(index); }
        [[nodiscard]] __device__ const ScalarType& operator[](size_t index) const { return *data_ptr(index); }
        /* Operations */
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return mat.getRow() * mat.getCol(); }
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(size_t index);
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(size_t index) const;
    };

    template<Matrix T>
    __host__ __device__ inline typename device_obj<ContinuousFlatten<T>>::PtrTy
    device_obj<ContinuousFlatten<T>>::data_ptr(size_t index) {
        const size_t major = index / mat.getMaxMinor();
        const size_t minor = index % mat.getMaxMinor();
        const size_t row = MatrixOption::rowFromMajorMinor<T>(major, minor);
        const size_t col = MatrixOption::colFromMajorMinor<T>(major, minor);
        return mat.data_ptr(row, col);
    }

    template<Matrix T>
    __host__ __device__ inline typename device_obj<ContinuousFlatten<T>>::ConstPtrTy
    device_obj<ContinuousFlatten<T>>::data_ptr(size_t index) const {
        return const_cast<This&>(*this).data_ptr(index);
    }
}

namespace Physica {
    template<Matrix T>
    class Traits<Core::device_obj<Core::ContinuousFlatten<T>>> : public Traits<Core::ContinuousFlatten<T>> {};
}
