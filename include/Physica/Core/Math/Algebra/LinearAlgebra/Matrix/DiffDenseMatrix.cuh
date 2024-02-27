/*
 * Copyright 2024 WeiBo He.
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

#include "DiffDenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/RValueMatrix.cuh"

namespace Physica::Core {
    namespace Internal {
        template<class T, int Option, unsigned int Order>
        class Traits<Core::device_obj<Differentiable<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>>
                : public Traits<DenseMatrix<T, Option>> {
        public:
            using ScalarType = device_obj<Differentiable<T, DiffMode::Reverse, Order>>;
        };
    }

    template<class PlainScalar, int Option, unsigned int Order>
    class device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>
            : public device_obj<RValueMatrix<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>>
            , public DenseMatrixDim<device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>, Dynamic, Dynamic, Dynamic, Dynamic> {
        static_assert(!PlainScalar::isDifferentiable, "[Error]: Nested Differentiable<> is not allowed");
        using PlainMatrix = DenseMatrix<PlainScalar, Option>;
        using host_obj = Differentiable<PlainMatrix, DiffMode::Reverse, Order>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using Dim = DenseMatrixDim<This, Dynamic, Dynamic, Dynamic, Dynamic>;
        using TracerType = device_obj<typename host_obj::TracerType>;
        using SegmentType = device_obj<typename host_obj::SegmentType>;
    public:
        using ScalarType = typename Base::ScalarType;
        static_assert(Utils::is_device_obj<ScalarType>::value, "Dbug");
        using DiffRecord = typename SegmentType::DiffRecord;
        using OperandArray = typename SegmentType::OperandArray;
    private:
        PlainStruct<SegmentType> traceSeg;
    public:
        device_obj() = default;
        device_obj(size_t row, size_t column, ExpressionType type);
        device_obj(const PlainMatrix& values);
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(const device_obj&) = default;
        device_obj& operator=(device_obj&&) noexcept = default;
        /* Operations */
        template<class RandomGenerator> inline void random_uniform(RandomGenerator& gen);
        template<class RandomGenerator> inline void random_normal(RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        inline void random_any(Distribution& dist, RandomGenerator& gen);

        [[nodiscard]] This copy() const;
        void swap(device_obj& obj) noexcept { std::swap(*this, obj); }
        /* Getters */
        using Dim::getRow;
        using Dim::getColumn;
        [[nodiscard]] __device__ inline size_t calcOffset(size_t row, size_t col) const;
        [[nodiscard]] __device__ inline ScalarType calc(size_t row, size_t col) const;
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return getTraceSegment().getLength(); }
        [[nodiscard]] __device__ inline DiffRecord& getRecord(size_t row, size_t col);
        [[nodiscard]] __device__ OperandArray& getOperands() { return getTraceSegment().getOperands(); }
        [[nodiscard]] __device__ inline PlainScalar& getValue(size_t row, size_t col);
        [[nodiscard]] __device__ inline const PlainScalar& getValue(size_t row, size_t col) const;
        [[nodiscard]] __device__ inline PlainScalar& getGrad(size_t row, size_t col);
        [[nodiscard]] __device__ inline const PlainScalar& getGrad(size_t row, size_t col) const;
        /* Static members */
        template<class RandomGenerator>
        [[nodiscard]] inline static This random_uniform(size_t row, size_t column, RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] inline static This random_normal(size_t row, size_t column, RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        [[nodiscard]] inline static This random_any(size_t row, size_t column, Distribution& dist, RandomGenerator& gen);
    private:
        /* Operations */
        using Dim::resize;
        /* Getters */
        [[nodiscard]] __host__ __device__ SegmentType& getTraceSegment() noexcept { return traceSeg.getDerived(); }
        [[nodiscard]] __host__ __device__ const SegmentType& getTraceSegment() const noexcept { return traceSeg.getDerived(); }
    };

    template<class PlainScalar, int Option, unsigned int Order>
    device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>::device_obj(size_t row, size_t column, ExpressionType type)
            : Dim(row, column)
            , traceSeg(asStruct(TracerType::getInstance().pushSegment(row * column, type))) {}

    template<class PlainScalar, int Option, unsigned int Order>
    device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>::device_obj(const PlainMatrix& values)
            : Dim(values.getRow(), values.getColumn())
            , traceSeg(asStruct(TracerType::getInstance().pushSegment(values.flatten()))) {}

    template<class PlainScalar, int Option, unsigned int Order>
    template<class RandomGenerator>
    inline void device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>::random_uniform(RandomGenerator& gen) {
        *this = random_uniform(getRow(), getColumn(), gen);
    }

    template<class PlainScalar, int Option, unsigned int Order>
    template<class RandomGenerator>
    inline void device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>::random_normal(RandomGenerator& gen) {
        *this = random_normal(getRow(), getColumn(), gen);
    }

    template<class PlainScalar, int Option, unsigned int Order>
    template<class Distribution, class RandomGenerator>
    inline void device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>::random_any(Distribution& dist, RandomGenerator& gen) {
        *this = random_any(getRow(), getColumn(), dist, gen);
    }

    template<class PlainScalar, int Option, unsigned int Order>
    device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>
    device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>::copy() const {
        This result{};
        const auto& newTrace = TracerType::getInstance().pushSegment(getTraceSegment().copy());
        result.traceSeg = asStruct(newTrace);
        return result;
    }

    template<class PlainScalar, int Option, unsigned int Order>
    __device__ inline size_t
    device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>::calcOffset(size_t row, size_t col) const {
        if constexpr (Base::isRowMatrix)
            return row * getColumn() + col;
        else
            return col * getRow() + row;
    }

    template<class PlainScalar, int Option, unsigned int Order>
    __device__ inline typename device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>::ScalarType
    device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>::calc(size_t row, size_t col) const {
        return getTraceSegment()[calcOffset(row, col)];
    }

    template<class PlainScalar, int Option, unsigned int Order>
    __device__ inline typename device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>::DiffRecord&
    device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>::getRecord(size_t row, size_t col) {
        return getTraceSegment().getRecords()[calcOffset(row, col)];
    }

    template<class PlainScalar, int Option, unsigned int Order>
    __device__ inline PlainScalar&
    device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>::getValue(size_t row, size_t col) {
        return getTraceSegment().getValues()[calcOffset(row, col)];
    }

    template<class PlainScalar, int Option, unsigned int Order>
    __device__ inline const PlainScalar&
    device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>::getValue(size_t row, size_t col) const {
        return getTraceSegment().getValues()[calcOffset(row, col)];
    }

    template<class PlainScalar, int Option, unsigned int Order>
    __device__ inline PlainScalar&
    device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>::getGrad(size_t row, size_t col) {
        return getTraceSegment().getGrads()[calcOffset(row, col)];
    }

    template<class PlainScalar, int Option, unsigned int Order>
    __device__ inline const PlainScalar&
    device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>::getGrad(size_t row, size_t col) const {
        return getTraceSegment().getGrads()[calcOffset(row, col)];
    }

    template<class PlainScalar, int Option, unsigned int Order>
    template<class RandomGenerator>
    inline device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>
    device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>::random_uniform(
            size_t row, size_t column, RandomGenerator& gen) {
        return This(PlainMatrix::random_uniform(row, column, gen));
    }

    template<class PlainScalar, int Option, unsigned int Order>
    template<class RandomGenerator>
    inline device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>
    device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>::random_normal(
            size_t row, size_t column, RandomGenerator& gen) {
        return This(PlainMatrix::random_normal(row, column, gen));
    }

    template<class PlainScalar, int Option, unsigned int Order>
    template<class Distribution, class RandomGenerator>
    inline device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>
    device_obj<Differentiable<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>::random_any(
            size_t row, size_t column, Distribution& dist, RandomGenerator& gen) {
        return This(PlainMatrix::random_any(row, column, dist, gen));
    }
}

#include "MatrixProduct/DiffMatrixProduct.cuh"
