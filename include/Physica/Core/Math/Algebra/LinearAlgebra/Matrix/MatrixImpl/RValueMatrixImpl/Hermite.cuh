/*
 * Copyright 2025 Weibo He.
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

#include "Hermite.h"

namespace Physica {
    template<class T>
    struct remove_hermite<device_obj<Hermite<T>>> {
        using Type = device_obj<T>;
    };

    template<Matrix M>
    class device_obj<Hermite<M>> : public device_obj<RValueMatrix<Hermite<M>>> {
        using host_obj = Hermite<M>;
        using Base = device_obj<RValueMatrix<host_obj>>;
    public:
        using typename Base::ScalarType;
    private:
        const device_obj<M>& matrix;
    public:
        explicit device_obj(const device_obj<M>& matrix_);
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const { return matrix.calc(col, row).conjugate(); }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return matrix.getCol(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return matrix.getRow(); }
    };

    template<Matrix M>
    device_obj<Hermite<M>>::device_obj(const device_obj<M>& matrix_) : matrix(matrix_) {
        static_assert(M::isComplex, "[Error]: Do not call hermite on real matrix");
    }

    template<Vector V>
    class device_obj<Hermite<V>> : public device_obj<RValueMatrix<Hermite<V>>> {
        using host_obj = Hermite<V>;
        using Base = device_obj<RValueMatrix<host_obj>>;
    public:
        using typename Base::ScalarType;
    private:
        const device_obj<V>& vec;
    public:
        explicit device_obj(const device_obj<V>& vec_);
        /* Operations */
        __device__ void assign(Matrix auto& target) const;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return vec.calc(col).conjugate(); }
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return vec.getLength(); }
    };

    template<Vector V>
    device_obj<Hermite<V>>::device_obj(const device_obj<V>& vec_) : vec(vec_) {
        static_assert(V::isComplex, "[Error]: Do not call hermite on real vector");
    }

    template<Vector V>
    __device__ void device_obj<Hermite<V>>::assign(Matrix auto& target) const {
        for (size_t i = 0; i < vec.getLength(); ++i)
            target.refFromMajorMinor(0, i) = calc(target.rowFromMajorMinor(0, i), target.colFromMajorMinor(0, i));
    }
}

namespace Physica {
    template <Matrix M>
    class Traits<device_obj<Hermite<M>>> : public Traits<Hermite<M>> {};

    template <Vector V>
    class Traits<device_obj<Hermite<V>>> : public Traits<Hermite<V>> {};
}

