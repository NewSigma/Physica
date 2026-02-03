/*
 * Copyright 2021-2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

namespace Physica {
    template<class T,
             int Option = MatrixOption::Col,
             size_t Row = Dynamic,
             size_t Col = Dynamic,
             class Allocator = HostAllocator<T>>
    class Array2D {
        using This = Array2D<T, Option, Row, Col, Allocator>;
        constexpr static bool isDynamicArray = Row == Dynamic && Col == Dynamic;
        constexpr static bool isColMajor = (Option & MatrixOption::Col) != 0;
        constexpr static size_t MaxMajor = isColMajor ? Col : Row;
        constexpr static size_t MaxMinor = isColMajor ? Row : Col;
        constexpr static int TransOption = isColMajor ? MatrixOption::Row : MatrixOption::Col;
        static_assert(Option != MatrixOption::BothMajor, "[Error]: Invalid Option");

        template<class U>
        struct Helper {
            using ArrayType = Array<T, Row * Col>;
        };

        template<Scalar U>
        struct Helper<U> {
            using ArrayType = DenseVector<T, Row * Col>;
        };
    private:
        using ArrayType = Helper<T>::ArrayType;
        using IndexType = std::conditional<isDynamicArray, size_t, PlainStruct<void>>::type;

        ArrayType arr;
        [[no_unique_address]] IndexType r = 0;
    public:
        Array2D() = default;
        explicit Array2D(size_t order);
        Array2D(size_t row, size_t col, auto&&... args);
        Array2D(std::initializer_list<T> list);
        Array2D(std::initializer_list<ArrayType> list);
        Array2D(const This&) = default;
        Array2D(This&&) noexcept = default;
        ~Array2D() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] T& operator[](size_t r, size_t c);
        [[nodiscard]] const T& operator[](size_t r, size_t c) const;
        /* Operations */
        void resize(size_t row, size_t col, auto&&... args);
        void resize(size_t order);
        [[nodiscard]] auto toDevice() const;
        [[nodiscard]] auto toDeviceAsync() const;
        void toDevice(device_obj<This>& obj) const;
        void toDeviceAsync(device_obj<This>& obj) const;

        [[nodiscard]] auto transpose() const noexcept;

        void zeros() noexcept;
        void read(const T* __restrict p) noexcept;
        void swap(This& __restrict obj) noexcept;
        void swap_row(size_t r1, size_t r2) noexcept;
        void swap_col(size_t c1, size_t c2) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ auto&& asArray(this auto&& self) noexcept { return self.arr; }
        [[nodiscard, gnu::returns_nonnull]] __host__ __device__ auto* data(this auto&&) noexcept;
        [[nodiscard, gnu::returns_nonnull]] __host__ __device__ auto* data_ptr(this auto&&, size_t row, size_t col) noexcept;
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept;
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept;
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept;
        [[nodiscard]] __host__ __device__ size_t getMaxMajor() const noexcept;
        [[nodiscard]] __host__ __device__ size_t getMaxMinor() const noexcept;
        [[nodiscard]] __host__ __device__ bool empty() const noexcept;
        /* Static members */
        [[nodiscard]] static This read(size_t row, size_t col, const T* __restrict p) noexcept;
    private:
        Array2D(ArrayType arr_, IndexType r_);

        [[nodiscard]] __host__ __device__ size_t toIndex1D(size_t r, size_t c) const noexcept;
        /* Friends */
        friend class Array2D<T, TransOption, Col, Row, Allocator>;
        friend class device_obj<This>;
    };
}

namespace Physica {
    template<class T, int Op, size_t Row, size_t Col, class Allocator>
    class Traits<Array2D<T, Op, Row, Col, Allocator>> {
    public:
        constexpr static int Option = Op;
    };
}

#include "ArrayImpl/Array2DImpl.h"
