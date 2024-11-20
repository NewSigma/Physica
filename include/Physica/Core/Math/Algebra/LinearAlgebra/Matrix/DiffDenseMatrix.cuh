/*
 * Copyright 2024 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/RValueMatrix.cuh"
#include "DiffDenseMatrix.h"

namespace Physica::Core {
    template<Scalar T, int Option, int Order>
    class device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>
            : public device_obj<RValueMatrix<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>>
            , public DenseMatrixDim<device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>, Dynamic, Dynamic> {
        static_assert(!T::isDifferentiable, "[Error]: Nested Diff<> is not allowed");
        using PlainMatrix = DenseMatrix<T, Option>;
        using host_obj = Diff<PlainMatrix, DiffMode::Reverse, Order>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using Dim = DenseMatrixDim<This, Dynamic, Dynamic>;
        using TracerType = device_obj<typename host_obj::TracerType>;
        using SegmentType = device_obj<typename host_obj::SegmentType>;
    public:
        using ScalarType = typename Base::ScalarType;
        static_assert(is_device_obj<ScalarType>::value, "Invalid ScalarType");
        using DiffRecord = typename SegmentType::DiffRecord;
        using OperandArray = typename SegmentType::OperandArray;
    private:
        PlainStruct<SegmentType> traceSeg;
    public:
        device_obj() = default;
        device_obj(size_t row, size_t col, ExprType type);
        device_obj(const PlainMatrix& values);
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(const device_obj&) = default;
        device_obj& operator=(device_obj&&) noexcept = default;
        /* Operations */
        [[nodiscard]] __device__ inline size_t calcOffset(size_t row, size_t col) const;
        template<Side Owner = GetSide()>
        [[nodiscard]] __device__ inline ScalarType calc(size_t row, size_t col) const;

        template<class RandomType> inline void random_uniform(RandomType& gen);
        template<class RandomType> inline void random_normal(RandomType& gen);
        template<class Distribution, class RandomType>
        inline void random_any(Distribution& dist, RandomType& gen);

        [[nodiscard]] This copy() const;
        void swap(device_obj& obj) noexcept { std::swap(*this, obj); }
        /* Getters */
        using Dim::getRow;
        using Dim::getCol;
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return getTraceSegment().getLength(); }
        [[nodiscard]] __device__ inline DiffRecord& getRecord(size_t row, size_t col);
        [[nodiscard]] __device__ OperandArray& getOperands() { return getTraceSegment().getOperands(); }
        [[nodiscard]] __device__ inline T& getValue(size_t row, size_t col);
        [[nodiscard]] __device__ inline const T& getValue(size_t row, size_t col) const;
        [[nodiscard]] __device__ inline T& getGrad(size_t row, size_t col);
        [[nodiscard]] __device__ inline const T& getGrad(size_t row, size_t col) const;
        /* Static members */
        template<class RandomType>
        [[nodiscard]] inline static This random_uniform(size_t row, size_t col, RandomType& gen);
        template<class RandomType>
        [[nodiscard]] inline static This random_normal(size_t row, size_t col, RandomType& gen);
        template<class Distribution, class RandomType>
        [[nodiscard]] inline static This random_any(size_t row, size_t col, Distribution& dist, RandomType& gen);
    private:
        /* Operations */
        using Dim::resize;
        /* Getters */
        [[nodiscard]] __host__ __device__ SegmentType& getTraceSegment() noexcept { return traceSeg.getDerived(); }
        [[nodiscard]] __host__ __device__ const SegmentType& getTraceSegment() const noexcept { return traceSeg.getDerived(); }
    };

    template<Scalar T, int Option, int Order>
    device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>::device_obj(size_t row, size_t col, ExprType type)
            : Dim(row, col)
            , traceSeg(asStruct(TracerType::getInstance().pushSegment(row * col, type))) {}

    template<Scalar T, int Option, int Order>
    device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>::device_obj(const PlainMatrix& values)
            : Dim(values.getRow(), values.getCol())
            , traceSeg(asStruct(TracerType::getInstance().pushSegment(values.flatten()))) {}

    template<Scalar T, int Option, int Order>
    __device__ inline size_t
    device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>::calcOffset(size_t row, size_t col) const {
        if constexpr (MatrixOption::isRowMatrix<This>())
            return row * getCol() + col;
        else
            return col * getRow() + row;
    }

    template<Scalar T, int Option, int Order>
    template<Side Owner>
    __device__ inline typename device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>::ScalarType
    device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>::calc(size_t row, size_t col) const {
        return getTraceSegment()[calcOffset(row, col)];
    }

    template<Scalar T, int Option, int Order>
    template<class RandomType>
    inline void device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>::random_uniform(RandomType& gen) {
        *this = random_uniform(getRow(), getCol(), gen);
    }

    template<Scalar T, int Option, int Order>
    template<class RandomType>
    inline void device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>::random_normal(RandomType& gen) {
        *this = random_normal(getRow(), getCol(), gen);
    }

    template<Scalar T, int Option, int Order>
    template<class Distribution, class RandomType>
    inline void device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>::random_any(Distribution& dist, RandomType& gen) {
        *this = random_any(getRow(), getCol(), dist, gen);
    }

    template<Scalar T, int Option, int Order>
    device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>
    device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>::copy() const {
        This result{};
        const auto& newTrace = TracerType::getInstance().pushSegment(getTraceSegment().copy());
        result.traceSeg = asStruct(newTrace);
        return result;
    }

    template<Scalar T, int Option, int Order>
    __device__ inline typename device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>::DiffRecord&
    device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>::getRecord(size_t row, size_t col) {
        return getTraceSegment().getRecords()[calcOffset(row, col)];
    }

    template<Scalar T, int Option, int Order>
    __device__ inline T&
    device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>::getValue(size_t row, size_t col) {
        return getTraceSegment().getValues()[calcOffset(row, col)];
    }

    template<Scalar T, int Option, int Order>
    __device__ inline const T&
    device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>::getValue(size_t row, size_t col) const {
        return getTraceSegment().getValues()[calcOffset(row, col)];
    }

    template<Scalar T, int Option, int Order>
    __device__ inline T&
    device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>::getGrad(size_t row, size_t col) {
        return getTraceSegment().getGrads()[calcOffset(row, col)];
    }

    template<Scalar T, int Option, int Order>
    __device__ inline const T&
    device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>::getGrad(size_t row, size_t col) const {
        return getTraceSegment().getGrads()[calcOffset(row, col)];
    }

    template<Scalar T, int Option, int Order>
    template<class RandomType>
    inline device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>
    device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>::random_uniform(
            size_t row, size_t col, RandomType& gen) {
        return This(PlainMatrix::random_uniform(row, col, gen));
    }

    template<Scalar T, int Option, int Order>
    template<class RandomType>
    inline device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>
    device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>::random_normal(
            size_t row, size_t col, RandomType& gen) {
        return This(PlainMatrix::random_normal(row, col, gen));
    }

    template<Scalar T, int Option, int Order>
    template<class Distribution, class RandomType>
    inline device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>
    device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>::random_any(
            size_t row, size_t col, Distribution& dist, RandomType& gen) {
        return This(PlainMatrix::random_any(row, col, dist, gen));
    }
}

namespace Physica {
    template<Scalar T, int Option, int Order>
    class Traits<Core::device_obj<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>>
            : public Traits<DenseMatrix<T, Option>> {
    public:
        using ScalarType = Core::device_obj<Diff<T, DiffMode::Reverse, Order>>;
    };
}

#include "MatrixImpl/MatrixProduct/DiffMatrixProduct.cuh"
