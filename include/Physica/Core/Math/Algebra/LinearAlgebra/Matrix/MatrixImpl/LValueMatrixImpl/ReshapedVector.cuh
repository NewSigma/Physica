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
    protected:
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
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
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept;
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept;
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(size_t row, size_t col) noexcept;
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(size_t row, size_t col) const noexcept;
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
    __host__ __device__ inline auto device_obj<LValueReshapedVector<T, MatrixMajor, Row, Col>>::data_ptr(size_t row, size_t col) noexcept -> PtrTy {
        assert(row < getRow() && col < getCol());
        if constexpr (MatrixOption::isColMatrix<This>())
            return v.data_ptr(col * getRow() + row);
        else
            return v.data_ptr(row * getCol() + col);
    }

    template<Vector T, int MatrixMajor, size_t Row, size_t Col>
    __host__ __device__ inline auto device_obj<LValueReshapedVector<T, MatrixMajor, Row, Col>>::data_ptr(size_t row, size_t col) const noexcept -> ConstPtrTy {
        return const_cast<This&>(*this).data_ptr(row, col);
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
