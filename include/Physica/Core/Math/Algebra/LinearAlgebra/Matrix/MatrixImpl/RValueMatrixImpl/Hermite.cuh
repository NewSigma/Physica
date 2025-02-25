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
    template<Matrix T>
    class device_obj<Hermite<T>> : public device_obj<RValueMatrix<Hermite<T>>> {
        using host_obj = Hermite<T>;
        using Base = device_obj<RValueMatrix<host_obj>>;
    public:
        using typename Base::ScalarType;
    private:
        const device_obj<T>& matrix;
    public:
        device_obj(const device_obj<T>& matrix_) : matrix(matrix_) {}
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const { return matrix.calc(col, row).conjugate(); }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return matrix.getCol(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return matrix.getRow(); }
    };

    template<Vector T>
    class device_obj<HermiteVector<T>> : public device_obj<RValueMatrix<HermiteVector<T>>> {
        using host_obj = HermiteVector<T>;
        using Base = device_obj<RValueMatrix<host_obj>>;
    public:
        using typename Base::ScalarType;
    private:
        const device_obj<T>& vec;
    public:
        explicit device_obj(const device_obj<T>& vec_) : vec(vec_) {}
        /* Operations */
        template<Matrix M>
        __device__ void assign(M& target) const;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return vec.calc(col).conjugate(); }
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return vec.getLength(); }
    };

    template<Vector T>
    template<Matrix M>
    __device__ void device_obj<HermiteVector<T>>::assign(M& target) const {
        for (size_t i = 0; i < vec.getLength(); ++i)
            target.refFromMajorMinor(0, i) = calc(M::rowFromMajorMinor(0, i), M::colFromMajorMinor(0, i));
    }
}

namespace Physica {
    template <Matrix T>
    class Traits<device_obj<Hermite<T>>> : public Traits<Hermite<T>> {};

    template <Vector T>
    class Traits<device_obj<HermiteVector<T>>> : public Traits<HermiteVector<T>> {};
}

