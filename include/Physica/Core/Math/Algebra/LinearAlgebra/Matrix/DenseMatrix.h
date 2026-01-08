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
#include "MatrixImpl/ContinuousMatrix.h"

namespace Physica {
    /**
     * \class DenseMatrix
     * A matrix can be either fixed matrix, which have its max size defined,
     * or dynamic matrix, whose size is dynamically changed.
     *
     * \tparam Option
     * Option is combinations of \enum MatrixOption
     */
    template<Scalar T,
             int Option = MatrixOption::Col,
             size_t Row = Dynamic,
             size_t Col = Dynamic,
             class Allocator = HostAllocator<T>>
    class DenseMatrix : public ContinuousMatrix<DenseMatrix<T, Option, Row, Col, Allocator>>
                      , public CRCoro<DenseMatrix<T, Option, Row, Col, Allocator>>
                      , public Array2D<T, Option, Row, Col, Allocator> {
        using This = DenseMatrix<T, Option, Row, Col, Allocator>;
        using Base = ContinuousMatrix<This>;
        using Storage = Array2D<T, Option, Row, Col, Allocator>;
        using Base::isReverseDiff;
    public:
        using typename Base::ScalarType;
        using device_obj_type = device_obj<This>;
        using ColMatrix = DenseMatrix<T, MatrixOption::Col, Row, Col>;
        using RowMatrix = DenseMatrix<T, MatrixOption::Row, Row, Col>;
        using VectorIniter = DenseVector<T, MatrixOption::isColMatrix<This>() ? Row : Col>;
        template<Scalar U>
        using rebind_scalar = DenseMatrix<U, Option, Row, Col, Allocator>;
    public:
        DenseMatrix() = default;
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
        using Base::operator[];
        /* Operations */
        size_t completePivoting(size_t col);
        size_t partialPivoting(size_t col);

        using Storage::resize;
        void resize(const Matrix auto& m, auto&&... args);
        [[nodiscard]] auto toDevice() const;
        [[nodiscard]] auto toDeviceAsync() const;
        using Base::toDevice;
        using Base::toDeviceAsync;

        using Base::transpose;
        [[nodiscard]] auto&& flatten(this auto&&) noexcept;

        using Base::random_any;
        using Base::random_normal;
        using Base::random_uniform;
        using Storage::zeros;
        using Storage::swap;

        using Base::read;
        /* Getters */
        using Storage::data;
        using Storage::data_ptr;
        using Storage::getCol;
        using Storage::getRow;
        using Storage::getSize;
        using Storage::getMaxMajor;
        using Storage::getMaxMinor;
        /* Static members */
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
    private:
        DenseMatrix(Storage storage) : Storage(std::move(storage)) {}
        friend class device_obj<This>;
    };

    // Just give me a matrix. This is what you'll get.
    template<Scalar T> using MatrixND = DenseMatrix<T>;

    template<Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    void swap(DenseMatrix<T, Option, Row, Col, Allocator>& __restrict m1,
              DenseMatrix<T, Option, Row, Col, Allocator>& __restrict m2) noexcept {
        m1.swap(m2);
    }
}

namespace Physica {
    template<Scalar T, int Op, size_t Row, size_t Col, class Allocator>
    class Traits<DenseMatrix<T, Op, Row, Col, Allocator>> {
        static_assert(!Diffable<T>, "[Error]: Use diffable matrix instead");
    public:
        using ScalarType = T;
        constexpr static int Option = Op;
        constexpr static size_t RowAtCompile = Row;
        constexpr static size_t ColAtCompile = Col;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
        using AllocatorType = Allocator;

        constexpr static bool FastAssign = false;
    };
}

namespace std {
    template<Physica::Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    struct formatter<Physica::DenseMatrix<T, Option, Row, Col, Allocator>, char> {
        constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }
        static auto format(const Physica::DenseMatrix<T, Option, Row, Col, Allocator>& obj, auto& ctx) {
            auto f = obj.format();
            return formatter<decltype(f), char>::format(f, ctx);
        }
    };
}

#include "MatrixImpl/DenseMatrixImpl.h"
