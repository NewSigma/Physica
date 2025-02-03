/*
 * Copyright 2023-2025 Weibo He.
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

#include "../RValueMatrix.cuh"

namespace Physica::Core {
    template<Matrix T>
    class device_obj<Transpose<T>> : public device_obj<RValueMatrix<Transpose<T>>> {
        using host_obj = Transpose<T>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;

        const device_obj<T>& mat;
    public:
        using typename Base::ScalarType;
    public:
        __host__ __device__ device_obj(const device_obj<T>& mat_) : mat(mat_) {}
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const { return mat.calc(col, row); }
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return mat.template getCol<Owner>(); }
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return mat.template getRow<Owner>(); }
    };

    template<Vector T>
    class device_obj<TransposeVector<T>> : public device_obj<RValueMatrix<TransposeVector<T>>> {
        using host_obj = TransposeVector<T>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;

        const device_obj<T>& vec;
    public:
        using typename Base::ScalarType;
    public:
        __host__ __device__ explicit device_obj(const device_obj<T>& vec_) : vec(vec_) {}
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return vec.calc(col); }
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return 1; }
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return vec.getLength(); }
    };
}

namespace Physica {
    template<Matrix T>
    class Traits<device_obj<Transpose<T>>> : public Traits<Transpose<T>> {};

    template<Vector T>
    class Traits<device_obj<TransposeVector<T>>> : public Traits<TransposeVector<T>> {};
}
