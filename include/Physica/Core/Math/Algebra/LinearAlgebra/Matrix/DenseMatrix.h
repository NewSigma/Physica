/*
 * Copyright 2020-2024 Weibo He.
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

#include "DenseMatrixImpl/DenseMatrixStorage.h"
#include "MatrixImpl/ContinuousMatrix.h"

namespace Physica::Core {
    /**
     * \class DenseMatrix
     * A matrix can be either fixed matrix, which have its max size defined,
     * or dynamic matrix, whose size is dynamically changed.
     *
     * \tparam Option
     * Option is combinations of \enum MatrixOption
     */
    template<Scalar T,
             int Option = MatrixOption::Col | MatrixOption::Vector,
             size_t Row = Dynamic,
             size_t Col = Dynamic,
             class Allocator = HostAllocator<T>>
    class DenseMatrix : public ContinuousMatrix<DenseMatrix<T, Option, Row, Col, Allocator>>
                      , public CRCoro<DenseMatrix<T, Option, Row, Col, Allocator>>
                      , public DenseMatrixStorage<T, Option, Row, Col, Allocator> {
        using This = DenseMatrix<T, Option, Row, Col, Allocator>;
        using Base = ContinuousMatrix<This>;
        using Storage = DenseMatrixStorage<T, Option, Row, Col, Allocator>;
        using Base::isReverseDiff;
    public:
        using typename Base::ValueType;
        using typename Base::ScalarType;
        using device_obj_type = device_obj<This>;
        using initializer_list = std::initializer_list<typename Storage::InitializerType>;
        using ColMatrix = DenseMatrix<T, MatrixOption::getStorage<DenseMatrix>() | MatrixOption::Col, Row, Col>;
        using RowMatrix = DenseMatrix<T, MatrixOption::getStorage<DenseMatrix>() | MatrixOption::Row, Row, Col>;
        using RealMatrix = DenseMatrix<typename T::RealType, Option, Row, Col>;
    public:
        DenseMatrix() = default;
        DenseMatrix(size_t row, size_t col);
        DenseMatrix(size_t row, size_t col, T value);
        DenseMatrix(initializer_list list);
        template<Matrix M>
        DenseMatrix(const M& mat);
        template<Vector V>
        DenseMatrix(const V& mat);
        DenseMatrix(const This&) = default;
        DenseMatrix(This&&) noexcept = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        using Base::operator=;
        using Base::operator();
        /* Operations */
        size_t completePivoting(size_t col);
        size_t partialPivoting(size_t col);

        using Storage::resize;
        [[nodiscard]] inline auto toDevice() const;
        [[nodiscard]] inline auto toDeviceAsync() const;
        void toDevice(device_obj<This>& obj) const;
        void toDeviceAsync(device_obj<This>& obj) const;
        using Storage::swap;

        using Base::random_any;
        using Base::random_normal;
        using Base::random_uniform;
        /* Getters */
        using Base::data;
        using Storage::data_ptr;
        using Storage::getCol;
        using Storage::getRow;
        /* Static members */
        [[nodiscard]] static DenseMatrix zeros(size_t rank) { return DenseMatrix(rank, rank, T(0)); }
        [[nodiscard]] static DenseMatrix zeros(size_t row, size_t col) { return DenseMatrix(row, col, T(0)); }
        [[nodiscard]] static DenseMatrix unitMatrix(size_t order);
        template<RandomGenerator R>
        [[nodiscard]] static DenseMatrix random_uniform(size_t order) { return random_uniform<R>(order, order); }
        template<RandomGenerator R>
        [[nodiscard]] inline static auto random_uniform(size_t row, size_t col);
        template<RandomGenerator R>
        [[nodiscard]] inline static auto random_normal(size_t row, size_t col);
        template<class Distribution, RandomGenerator R>
        [[nodiscard]] inline static auto random_any(size_t row, size_t col, Distribution& dist);
        template<Vector V>
        [[nodiscard]] static std::pair<DenseMatrix, DenseMatrix> meshgrid(const V& vecInCols, const V& vecInRows);
    private:
        DenseMatrix(Storage storage) : Storage(std::move(storage)) {}
        friend class device_obj<This>;
    };
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
    };
}

namespace std {
    template<Physica::Core::Scalar T, int Option, size_t Row, size_t Col, class Allocator>
    inline void swap(
            Physica::Core::DenseMatrix<T, Option, Row, Col, Allocator>& __restrict m1,
            Physica::Core::DenseMatrix<T, Option, Row, Col, Allocator>& __restrict m2) noexcept {
        m1.swap(m2);
    }
}

#include "DenseMatrixImpl/DenseMatrixImpl.h"
