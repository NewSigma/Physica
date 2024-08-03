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

#include <Physica/Core/MultiPrecision/Differentiable.cuh>
#include "DiffVector.h"
#include "Vector.cuh"

namespace Physica::Core {
    template<class PlainScalar, unsigned Order>
    class device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse, Order>>
            : public device_obj<RValueVector<Differentiable<Vector<PlainScalar>, DiffMode::Reverse, Order>>> {
        static_assert(!PlainScalar::isDifferentiable, "[Error]: Nested Differentiable<> is not allowed");
        using PlainVector = Vector<PlainScalar>;
        using host_obj = Differentiable<PlainVector, DiffMode::Reverse, Order>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
        using TracerType = device_obj<typename host_obj::TracerType>;
        using SegmentType = device_obj<typename host_obj::SegmentType>;
        using RecordArray = typename SegmentType::RecordArray;
        using ValueVector = typename SegmentType::ValueVector;
        using GradVector = typename SegmentType::GradVector;
    public:
        using ScalarType = typename Base::ScalarType;
        using DiffRecord = typename SegmentType::DiffRecord;
        using OperandArray = typename SegmentType::OperandArray;
        using Base::MaxThreadPerBlock;
    private:
        PlainStruct<SegmentType> traceSeg;
    public:
        device_obj() = default;
        device_obj(size_t length, ExpressionType type);
        device_obj(const PlainVector& values);
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

        template<bool ComputeMax>
        __device__ void minmaxKernelImpl(SegmentType& result) const;
        __device__ void sumKernelImpl(SegmentType& result) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ inline ScalarType calc(size_t index) const;
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return getTraceSegment().getLength(); }
        [[nodiscard]] __host__ __device__ inline PlainScalar* value_ptr(size_t index) const noexcept;
        [[nodiscard]] __host__ __device__ inline PlainScalar* grad_ptr(size_t index) const noexcept;
        [[nodiscard]] __device__ inline DiffRecord& getRecord(size_t index);
        [[nodiscard]] __device__ OperandArray& getOperands() noexcept { return getTraceSegment().getOperands(); }
        [[nodiscard]] __device__ inline PlainScalar& getValue(size_t index);
        [[nodiscard]] __device__ inline const PlainScalar& getValue(size_t index) const;
        [[nodiscard]] __device__ inline PlainScalar& getGrad(size_t index);
        [[nodiscard]] __device__ inline const PlainScalar& getGrad(size_t index) const;
        [[nodiscard]] const ValueVector& getValues() const noexcept { return getTraceSegment().getValues(); }
        [[nodiscard]] const GradVector& getGrads() const noexcept { return getTraceSegment().getGrads(); }
        [[nodiscard]] ScalarType max() const;
        [[nodiscard]] ScalarType min() const;
        [[nodiscard]] ScalarType sum() const;
        /* Static members */
        template<class RandomGenerator>
        [[nodiscard]] inline static This random_uniform(size_t len, RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] inline static This random_normal(size_t len, RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        [[nodiscard]] inline static This random_any(size_t len, Distribution& dist, RandomGenerator& gen);
    private:
        [[nodiscard]] __host__ __device__ SegmentType& getTraceSegment() noexcept { return traceSeg.getDerived(); }
        [[nodiscard]] __host__ __device__ const SegmentType& getTraceSegment() const noexcept { return traceSeg.getDerived(); }
    };
}

namespace Physica {
    using namespace Core;

    template<class T, unsigned int Order>
    class Traits<Core::device_obj<Differentiable<Vector<T>, DiffMode::Reverse, Order>>> : public Traits<Vector<T>> {
    public:
        using ScalarType = Core::device_obj<Differentiable<T, DiffMode::Reverse, Order>>;
    };
}

#include "DiffVectorImpl/DiffVectorImpl.cuh"
#include "DiffVectorImpl/DiffVectorExpr.cuh"
