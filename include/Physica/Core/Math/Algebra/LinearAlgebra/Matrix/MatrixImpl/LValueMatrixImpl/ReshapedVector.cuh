/*
 * Copyright 2025-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/LValueVector.cuh"
#include "../LValueMatrix.cuh"

namespace Physica {
    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    class device_obj<LValueReshapedVector<V, MatrixMajor, Row, Col>>
            : public device_obj<LValueMatrix<LValueReshapedVector<V, MatrixMajor, Row, Col>>> {
        using host_obj = LValueReshapedVector<V, MatrixMajor, Row, Col>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueMatrix<host_obj>>;
        using Ref = add_device_obj<V>::type;
    public:
        using typename Base::ScalarType;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<V>>> v;
        size_t r;
        size_t c;
    public:
    __host__ __device__ device_obj(Ref v_, size_t r_, size_t c_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        __host__ __device__ void resize(size_t row, size_t col);
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept;
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept;
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, size_t row, size_t col) noexcept;
        [[nodiscard]] ScalarType sum() const { return v.getDerived().sum(); }
    };

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    __host__ __device__ device_obj<LValueReshapedVector<V, MatrixMajor, Row, Col>>::device_obj(Ref v_, size_t r_, size_t c_)
            : v(asStruct(v_)), r(r_), c(c_) {
        assert(r == Row || Row == Dynamic);
        assert(c == Col || Col == Dynamic);
        assert(r * c == v.getLength());
    }

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    __host__ __device__ void device_obj<LValueReshapedVector<V, MatrixMajor, Row, Col>>::resize(
            [[maybe_unused]] size_t row, [[maybe_unused]] size_t col) {
        assert(row == getRow());
        assert(col == getCol());
    }

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    __host__ __device__ size_t device_obj<LValueReshapedVector<V, MatrixMajor, Row, Col>>::getRow() const noexcept {
        if constexpr (Row != Dynamic)
            return Row;
        return r;
    }

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    __host__ __device__ size_t device_obj<LValueReshapedVector<V, MatrixMajor, Row, Col>>::getCol() const noexcept {
        if constexpr (Col != Dynamic)
            return Col;
        return c;
    }

    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    __host__ __device__ auto device_obj<LValueReshapedVector<V, MatrixMajor, Row, Col>>::data_ptr(this auto&& self, size_t row, size_t col) noexcept {
        assert(row < self.getRow() && col < self.getCol());
        if constexpr (MatrixMajor::isColMatrix<This>())
            return self.v.getDerived().data_ptr(col * self.getRow() + row);
        else
            return self.v.getDerived().data_ptr(row * self.getCol() + col);
    }
}

namespace Physica {
    template<Vector V, int MatrixMajor, size_t Row, size_t Col>
    class Traits<device_obj<LValueReshapedVector<V, MatrixMajor, Row, Col>>>
            : public Traits<RValueReshapedVector<V, MatrixMajor, Row, Col>> {};
}
