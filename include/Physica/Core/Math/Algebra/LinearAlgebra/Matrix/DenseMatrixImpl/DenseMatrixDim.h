/*
 * Copyright 2022-2024 Weibo He.
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

#include <cassert>
#include <Physica/CRTPBase.h>

namespace Physica::Core {
    template<class Derived, size_t Row, size_t Col>
    class DenseMatrixDim {
        using This = DenseMatrixDim<Derived, Row, Col>;
    public:
        DenseMatrixDim() = default;
        __host__ __device__ DenseMatrixDim(size_t row, size_t col) { resize(row, col); }
        DenseMatrixDim(const This&) = default;
        DenseMatrixDim(This&&) noexcept = default;
        ~DenseMatrixDim() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Getters */
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return Row; }
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ constexpr static size_t getCol() noexcept { return Col; }
        /* Operations */
        __host__ __device__ void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) {
            assert(row == Row && col == Col && "Cannot resize a fixed matrix");
        }
        /* Helper */
        __host__ __device__ void swap([[maybe_unused]] This& dim) noexcept { /* Do nothing */ }
    };

    template<class Derived, size_t Col>
    class DenseMatrixDim<Derived, Dynamic, Col> : public CRTPBase<DenseMatrixDim<Derived, Dynamic, Col>> {
        using This = DenseMatrixDim<Derived, Dynamic, Col>;
        using Base = CRTPBase<This>;
    public:
        DenseMatrixDim() = default;
        __host__ __device__ DenseMatrixDim(size_t row, size_t col) { resize(row, col); }
        DenseMatrixDim(const This&) = default;
        DenseMatrixDim(This&&) noexcept = default;
        ~DenseMatrixDim() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Getters */
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept {
            const size_t size = Base::getDerived().getSize();
            assert(size % getCol() == 0);
            return size / getCol();
        }
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ constexpr static size_t getCol() noexcept { return Col; }
        /* Operations */
        __host__ __device__ void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) {
            assert(col == Col && "Cannot resize a fixed matrix");
        }
        /* Helper */
        __host__ __device__ void swap([[maybe_unused]] This& dim) noexcept { /* Do nothing */ }
    };

    template<class Derived, size_t Row>
    class DenseMatrixDim<Derived, Row, Dynamic> : public CRTPBase<DenseMatrixDim<Derived, Row, Dynamic>> {
        using This = DenseMatrixDim<Derived, Row, Dynamic>;
        using Base = CRTPBase<This>;
    public:
        DenseMatrixDim() = default;
        __host__ __device__ DenseMatrixDim(size_t row, size_t col) { resize(row, col); }
        DenseMatrixDim(const This&) = default;
        DenseMatrixDim(This&&) noexcept = default;
        ~DenseMatrixDim() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Getters */
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return Row; }
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept {
            const size_t size = Base::getDerived().getSize();
            assert(size % getRow() == 0);
            return size / getRow();
        }
        /* Operations */
        __host__ __device__ void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) {
            assert(row == Row && "Cannot resize a fixed matrix");
        }
        /* Helper */
        __host__ __device__ void swap([[maybe_unused]] This& dim) noexcept { /* Do nothing */ }
    };

    template<class Derived>
    class DenseMatrixDim<Derived, Dynamic, Dynamic> : public CRTPBase<DenseMatrixDim<Derived, Dynamic, Dynamic>> {
        using This = DenseMatrixDim<Derived, Dynamic, Dynamic>;
        using Base = CRTPBase<This>;
    private:
        size_t r;
    public:
        __host__ __device__ DenseMatrixDim() : r(0) {}
        __host__ __device__ DenseMatrixDim(size_t row, size_t col) { resize(row, col); }
        DenseMatrixDim(const This&) = default;
        DenseMatrixDim(This&&) noexcept = default;
        ~DenseMatrixDim() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Getters */
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return r; }
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept {
            const size_t size = Base::getDerived().getSize();
            assert(r == 0 || size % getRow() == 0);
            return r == 0 ? 0 : size / getRow();
        }
        /* Operations */
        __host__ __device__ void resize(size_t row, [[maybe_unused]] size_t col) { r = row; }
        /* Helper */
        __host__ __device__ void swap(This& __restrict dim) noexcept {
            assert(this != &dim && "[Error]: Self swap is likely a bug");
        #ifdef __CUDA_ARCH__
            thrust::swap(r, dim.r);
        #else
            std::swap(r, dim.r);
        #endif
        }
    };
}

namespace Physica {
    template<class T, size_t Row, size_t Col>
    class Traits<Core::DenseMatrixDim<T, Row, Col>> {
    public:
        using Derived = T;
    };
}
