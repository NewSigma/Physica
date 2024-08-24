/*
 * Copyright 2022 Weibo He.
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
#include "Physica/Utils/Template/CRTPBase.h"

namespace Physica::Core {
    template<class Derived, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixDim {
    public:
        DenseMatrixDim() = default;
        __host__ __device__ DenseMatrixDim(size_t row, size_t column) { resize(row, column); }
        DenseMatrixDim(const DenseMatrixDim&) = default;
        DenseMatrixDim(DenseMatrixDim&&) noexcept = default;
        ~DenseMatrixDim() = default;
        /* Operators */
        DenseMatrixDim& operator=(DenseMatrixDim obj) noexcept { swap(obj); return *this; }
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return Row; }
        [[nodiscard]] __host__ __device__ constexpr static size_t getColumn() noexcept { return Column; }
        /* Operations */
        __host__ __device__ void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t column) {
            assert(row == Row && column == Column && "Cannot resize a fixed matrix");
        }
        /* Helper */
        __host__ __device__ void swap([[maybe_unused]] DenseMatrixDim& dim) noexcept { /* Do nothing */ }
    };

    template<class Derived, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixDim<Derived, Dynamic, Column, MaxRow, MaxColumn> : public Utils::CRTPBase<Derived, 2> {
        using Base = Utils::CRTPBase<Derived, 2>;
    public:
        DenseMatrixDim() = default;
        __host__ __device__ DenseMatrixDim(size_t row, size_t column) { resize(row, column); }
        DenseMatrixDim(const DenseMatrixDim&) = default;
        DenseMatrixDim(DenseMatrixDim&&) noexcept = default;
        ~DenseMatrixDim() = default;
        /* Operators */
        DenseMatrixDim& operator=(DenseMatrixDim obj) noexcept { swap(obj); return *this; }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept {
            const size_t size = Base::getDerived().getSize();
            assert(size % getColumn() == 0);
            return size / getColumn();
        }
        [[nodiscard]] __host__ __device__ constexpr static size_t getColumn() noexcept { return Column; }
        /* Operations */
        __host__ __device__ void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t column) {
            assert(column == Column && "Cannot resize a fixed matrix");
        }
        /* Helper */
        __host__ __device__ void swap([[maybe_unused]] DenseMatrixDim& dim) noexcept { /* Do nothing */ }
    };

    template<class Derived, size_t Row, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixDim<Derived, Row, Dynamic, MaxRow, MaxColumn> : public Utils::CRTPBase<Derived, 2> {
        using Base = Utils::CRTPBase<Derived, 2>;
    public:
        DenseMatrixDim() = default;
        __host__ __device__ DenseMatrixDim(size_t row, size_t column) { resize(row, column); }
        DenseMatrixDim(const DenseMatrixDim&) = default;
        DenseMatrixDim(DenseMatrixDim&&) noexcept = default;
        ~DenseMatrixDim() = default;
        /* Operators */
        DenseMatrixDim& operator=(DenseMatrixDim obj) noexcept { swap(obj); return *this; }
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return Row; }
        [[nodiscard]] __host__ __device__ size_t getColumn() const noexcept {
            const size_t size = Base::getDerived().getSize();
            assert(size % getRow() == 0);
            return size / getRow();
        }
        /* Operations */
        __host__ __device__ void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t column) {
            assert(row == Row && "Cannot resize a fixed matrix");
        }
        /* Helper */
        __host__ __device__ void swap([[maybe_unused]] DenseMatrixDim& dim) noexcept { /* Do nothing */ }
    };

    template<class Derived, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixDim<Derived, Dynamic, Dynamic, MaxRow, MaxColumn> : public Utils::CRTPBase<Derived, 2> {
        using Base = Utils::CRTPBase<Derived, 2>;
    private:
        size_t r;
    public:
        __host__ __device__ DenseMatrixDim() : r(0) {}
        __host__ __device__ DenseMatrixDim(size_t row, size_t column) { resize(row, column); }
        DenseMatrixDim(const DenseMatrixDim&) = default;
        DenseMatrixDim(DenseMatrixDim&&) noexcept = default;
        ~DenseMatrixDim() = default;
        /* Operators */
        DenseMatrixDim& operator=(DenseMatrixDim obj) noexcept { swap(obj); return *this; }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return r; }
        [[nodiscard]] __host__ __device__ size_t getColumn() const noexcept {
            const size_t size = Base::getDerived().getSize();
            assert(r == 0 || size % getRow() == 0);
            return r == 0 ? 0 : size / getRow();
        }
        /* Operations */
        __host__ __device__ void resize(size_t row, [[maybe_unused]] size_t column) { r = row; }
        /* Helper */
        __host__ __device__ void swap(DenseMatrixDim& __restrict dim) noexcept {
            assert(this != &dim && "[Error]: Self swap is likely a bug");
        #ifdef __CUDA_ARCH__
            thrust::swap(r, dim.r);
        #else
            std::swap(r, dim.r);
        #endif
        }
    };
}
