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
    template<Vector T, int MatrixMajor, size_t Row, size_t Col>
    class device_obj<LValueReshapedVector<T, MatrixMajor, Row, Col>>
            : public device_obj<LValueMatrix<LValueReshapedVector<T, MatrixMajor, Row, Col>>> {
        using host_obj = LValueReshapedVector<T, MatrixMajor, Row, Col>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueMatrix<host_obj>>;
    public:
        using typename Base::ScalarType;
    private:
        device_obj<T>& v;
        size_t r;
        size_t c;
    public:
    __host__ __device__ device_obj(device_obj<T>& v_, size_t r_, size_t c_);
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
        [[nodiscard]] ScalarType sum() const { return v.sum(); }
    };

    template<Vector T, int MatrixMajor, size_t Row, size_t Col>
    __host__ __device__ device_obj<LValueReshapedVector<T, MatrixMajor, Row, Col>>::device_obj(device_obj<T>& v_, size_t r_, size_t c_)
            : v(v_), r(r_), c(c_) {
        assert(r == Row || Row == Dynamic);
        assert(c == Col || Col == Dynamic);
        assert(r * c == v.getLength());
    }

    template<Vector T, int MatrixMajor, size_t Row, size_t Col>
    __host__ __device__ void device_obj<LValueReshapedVector<T, MatrixMajor, Row, Col>>::resize(
            [[maybe_unused]] size_t row, [[maybe_unused]] size_t col) {
        assert(row == getRow());
        assert(col == getCol());
    }

    template<Vector T, int MatrixMajor, size_t Row, size_t Col>
    __host__ __device__ size_t device_obj<LValueReshapedVector<T, MatrixMajor, Row, Col>>::getRow() const noexcept {
        if constexpr (Row != Dynamic)
            return Row;
        return r;
    }

    template<Vector T, int MatrixMajor, size_t Row, size_t Col>
    __host__ __device__ size_t device_obj<LValueReshapedVector<T, MatrixMajor, Row, Col>>::getCol() const noexcept {
        if constexpr (Col != Dynamic)
            return Col;
        return c;
    }

    template<Vector T, int MatrixMajor, size_t Row, size_t Col>
    __host__ __device__ auto device_obj<LValueReshapedVector<T, MatrixMajor, Row, Col>>::data_ptr(this auto&& self, size_t row, size_t col) noexcept {
        assert(row < self.getRow() && col < self.getCol());
        if constexpr (MatrixOption::isColMatrix<This>())
            return self.v.data_ptr(col * self.getRow() + row);
        else
            return self.v.data_ptr(row * self.getCol() + col);
    }

    template<class Derived>
    template<Matrix M>
    __host__ __device__ auto device_obj<LValueVector<Derived>>::reshape(const M& mat) noexcept {
        using ResultType = device_obj<LValueReshapedVector<Derived, MatrixOption::getMajor<M>(), M::RowAtCompile, M::ColAtCompile>>;
        return ResultType(Base::getDerived(), mat.getRow(), mat.getCol());
    }

    template<class Derived>
    template<Matrix M>
    __host__ __device__ const auto device_obj<LValueVector<Derived>>::reshape(const M& mat) const noexcept {
        using ResultType = device_obj<LValueReshapedVector<Derived, MatrixOption::getMajor<M>(), M::RowAtCompile, M::ColAtCompile>>;
        return ResultType(Base::getConstCastDerived(), mat.getRow(), mat.getCol());
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ auto device_obj<LValueVector<Derived>>::reshape_col(size_t row, size_t col) noexcept {
        return device_obj<LValueReshapedVector<Derived, MatrixOption::Col, Row, Col>>(Base::getDerived(), row, col);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ const auto device_obj<LValueVector<Derived>>::reshape_col(size_t row, size_t col) const noexcept {
        return device_obj<LValueReshapedVector<Derived, MatrixOption::Col, Row, Col>>(Base::getConstCastDerived(), row, col);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ auto device_obj<LValueVector<Derived>>::reshape_row(size_t row, size_t col) noexcept {
        return device_obj<LValueReshapedVector<Derived, MatrixOption::Row, Row, Col>>(Base::getDerived(), row, col);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    __host__ __device__ const auto device_obj<LValueVector<Derived>>::reshape_row(size_t row, size_t col) const noexcept {
        return device_obj<LValueReshapedVector<Derived, MatrixOption::Row, Row, Col>>(Base::getConstCastDerived(), row, col);
    }
}

namespace Physica {
    template<Vector T, int MatrixMajor, size_t Row, size_t Col>
    class Traits<device_obj<LValueReshapedVector<T, MatrixMajor, Row, Col>>>
            : public Traits<RValueReshapedVector<T, MatrixMajor, Row, Col>> {};
}
