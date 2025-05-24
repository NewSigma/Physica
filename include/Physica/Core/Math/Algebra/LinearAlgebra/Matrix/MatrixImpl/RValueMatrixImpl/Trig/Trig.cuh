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

#include "../RValueMatrixImpl.cuh"
#include "Trig.h"

namespace Physica {
    template<Matrix M>
    class device_obj<TrigUpper<M>> : public device_obj<RValueMatrix<TrigUpper<M>>> {
        using host_obj = TrigUpper<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;

        device_obj<M>& mat;
    public:
        using typename Base::T;
        using typename Base::Tv;
    public:
        __host__ __device__ device_obj(device_obj<M>& mat_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const;
        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const;

        [[nodiscard]] T lnAbsDet() const;
        /* Getters */
        [[nodiscard]] __host__ __device__ const auto& getExpr() const noexcept { return mat; }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return mat.getCol(); }
    };

    template<Matrix M>
    __host__ __device__ device_obj<TrigUpper<M>>::device_obj(device_obj<M>& mat_) : mat(mat_) {
        assert(mat.isSquare());
    }

    template<Matrix M>
    __device__ auto device_obj<TrigUpper<M>>::calc(size_t row, size_t col) const -> T {
        if (row > col)
            return T(0);
        return mat.calc(row, col);
    }

    template<Matrix M>
    __device__ auto device_obj<TrigUpper<M>>::calc_value(size_t row, size_t col) const -> Tv {
        if (row > col)
            return Tv(0);
        return mat.calc_value(row, col);
    }

    template<Matrix M>
    auto device_obj<TrigUpper<M>>::lnAbsDet() const -> T {
        assert(!Base::diag().prod().isZero() && "[Error]: Singular matrix");
        return ln(abs(Base::diag())).sum();
    }

    template<Matrix M>
    class device_obj<TrigLower<M>> : public device_obj<RValueMatrix<TrigLower<M>>> {
        using host_obj = TrigLower<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;

        device_obj<M>& mat;
    public:
        using typename Base::T;
        using typename Base::Tv;
    public:
        __host__ __device__ device_obj(device_obj<M>& mat_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const;
        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const;

        [[nodiscard]] T lnAbsDet() const;
        /* Getters */
        [[nodiscard]] __host__ __device__ const auto& getExpr() const noexcept { return mat; }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return mat.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return mat.getCol(); }
    };

    template<Matrix M>
    __host__ __device__ device_obj<TrigLower<M>>::device_obj(device_obj<M>& mat_) : mat(mat_) {
        assert(mat.isSquare());
    }

    template<Matrix M>
    __device__ auto device_obj<TrigLower<M>>::calc(size_t row, size_t col) const -> T {
        if (row < col)
            return T(0);
        return mat.calc(row, col);
    }

    template<Matrix M>
    __device__ auto device_obj<TrigLower<M>>::calc_value(size_t row, size_t col) const -> Tv {
        if (row < col)
            return Tv(0);
        return mat.calc_value(row, col);
    }

    template<Matrix M>
    auto device_obj<TrigLower<M>>::lnAbsDet() const -> T {
        assert(!Base::diag().prod().isZero() && "[Error]: Singular matrix");
        return ln(abs(Base::diags())).sum();
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<device_obj<TrigUpper<M>>> : public Traits<TrigUpper<M>> {};

    template<Matrix M>
    class Traits<device_obj<TrigLower<M>>> : public Traits<TrigLower<M>> {};
}

#include "TrigGEMM.cuh"
