/*
 * Copyright 2020-2026 Weibo He.
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

#include "Physica/Core/Utils/Container/Array2D.h"
#include "MatrixImpl/CompactMatrix.h"

namespace Physica {
    template<Scalar T,
             int Major = MatrixMajor::Col,
             size_t Row = Dynamic,
             size_t Col = Dynamic,
             class Allocator = HostAllocator<T>>
    class DenseMatrix : public CompactMatrix<DenseMatrix<T, Major, Row, Col, Allocator>>
                      , public CRCoro<DenseMatrix<T, Major, Row, Col, Allocator>> {
        using This = DenseMatrix<T, Major, Row, Col, Allocator>;
        using Base = CompactMatrix<This>;
        using Storage = Array2D<T, Major, Row, Col, Allocator>;
        using Base::isReverseDiff;
    public:
        using typename Base::ScalarType;
        using device_obj_type = device_obj<This>;
        using ColMatrix = DenseMatrix<T, MatrixMajor::Col, Row, Col>;
        using RowMatrix = DenseMatrix<T, MatrixMajor::Row, Row, Col>;
        using VectorIniter = DenseVector<T, (Major == MatrixMajor::Col) ? Row : Col>;
        template<Scalar U>
        using rebind_scalar = DenseMatrix<U, Major, Row, Col, Allocator>;
    private:
        Storage storage;
    public:
        DenseMatrix() = default;
        explicit DenseMatrix(Storage storage) noexcept;
        explicit DenseMatrix(size_t order);
        DenseMatrix(size_t row, size_t col);
        DenseMatrix(size_t row, size_t col, T value);
        DenseMatrix(std::initializer_list<T> list);
        DenseMatrix(std::initializer_list<VectorIniter> list);
        DenseMatrix(const Matrix auto& mat);
        DenseMatrix(const Vector auto& vec);
        DenseMatrix(const This&) = default;
        DenseMatrix(This&&) noexcept = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        using Base::operator=;
        /* Operations */
        void resize(size_t order);
        void resize(const Matrix auto& m, auto&&... args);
        void resize(size_t row, size_t col, auto&&... args);

        [[nodiscard]] auto toDevice() const;
        [[nodiscard]] auto toDeviceAsync() const;
        using Base::toDevice;
        using Base::toDeviceAsync;

        [[nodiscard]] auto&& flatten(this auto&&) noexcept;

        using Base::random_any;
        using Base::random_normal;
        using Base::random_uniform;

        using Base::read;
        void zeros() noexcept;
        void swap_row(size_t r1, size_t r2) noexcept;
        void swap_col(size_t c1, size_t c2) noexcept;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] auto data(this auto&&) noexcept;
        [[nodiscard]] size_t getRow() const noexcept { return storage.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return storage.getCol(); }
        [[nodiscard]] size_t getSize() const noexcept { return storage.getSize(); }
        [[nodiscard]] auto&& asArray(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getRowAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getColAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept;
        [[nodiscard]] static This zeros(size_t order) { return zeros(order, order); }
        [[nodiscard]] static This zeros(size_t row, size_t col);
        [[nodiscard]] static This identity(size_t order);
        template<RNG R>
        [[nodiscard]] static This random_uniform(size_t order) { return random_uniform<R>(order, order); }
        template<RNG R>
        [[nodiscard]] static This random_uniform(size_t row, size_t col);
        template<RNG R>
        [[nodiscard]] static This random_normal(size_t row, size_t col);
        template<RNG R>
        [[nodiscard]] static This random_any(size_t row, size_t col, auto& distribution);
        [[nodiscard]] static auto meshgrid(const Vector auto& vecX, const Vector auto& vecY) -> std::pair<This, This>;
        [[nodiscard]] static This read(size_t row, size_t col, const T* __restrict p) noexcept;
        /* Friends */
        friend class device_obj<This>;
    };

    // Just give me a matrix. This is what you'll get.
    template<Scalar T> using MatrixND = DenseMatrix<T>;

    template<Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    void swap(DenseMatrix<T, Major, Row, Col, Allocator>& __restrict m1,
              DenseMatrix<T, Major, Row, Col, Allocator>& __restrict m2) noexcept {
        m1.swap(m2);
    }
}

namespace Physica {
    template<Scalar T, int Major_, size_t Row, size_t Col, class Allocator>
    class Traits<DenseMatrix<T, Major_, Row, Col, Allocator>> {
        static_assert(!Diffable<T>, "[Error]: Use diffable matrix instead");
    public:
        using ScalarType = T;
        using AllocatorType = Allocator;
    };
}

namespace std {
    template<Physica::Scalar T, int Major, size_t Row, size_t Col, class Allocator>
    struct formatter<Physica::DenseMatrix<T, Major, Row, Col, Allocator>, char> {
        constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }
        static auto format(const Physica::DenseMatrix<T, Major, Row, Col, Allocator>& obj, auto& ctx) {
            auto f = obj.format();
            return formatter<decltype(f), char>::format(f, ctx);
        }
    };
}

#include "MatrixImpl/DenseMatrixImpl.h"
