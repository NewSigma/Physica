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

#include "Physica/Core/MultiPrecision/Differentiable.cuh"
#include "DiffVector.h"
#include "Vector.cuh"

namespace Physica::Core {
    namespace Internal {
        template<class T>
        class Traits<device_obj<Differentiable<Vector<T>, DiffMode::Reverse>>> : public Traits<Vector<T>> {
        public:
            using ScalarType = device_obj<Differentiable<T, DiffMode::Reverse>>;
        };
    }

    template<class PlainScalar>
    class device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>
            : public device_obj<RValueVector<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>> {
        static_assert(!PlainScalar::isDifferentiable, "[Error]: Nested Differentiable<> is not allowed");
        using PlainVector = Vector<PlainScalar>;
        using host_obj = Differentiable<PlainVector, DiffMode::Reverse>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
        using TracerType = device_obj<typename host_obj::TracerType>;
        using SegmentType = device_obj<typename host_obj::SegmentType>;
        using RecordArray = typename SegmentType::RecordArray;
        using OperandArray = typename SegmentType::OperandArray;
        using DeviceVector = typename SegmentType::DeviceVector;
    public:
        using ScalarType = typename Base::ScalarType;
        using DiffScalar = typename SegmentType::DiffScalar;
        using DiffRecord = typename SegmentType::DiffRecord;
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
        [[nodiscard]] __device__ RecordArray& getRecords() noexcept { return getTraceSegment().getRecords(); }
        [[nodiscard]] __device__ OperandArray& getOperands() noexcept { return getTraceSegment().getOperands(); }
        [[nodiscard]] __device__ DeviceVector& getValues() noexcept { return values; }
        [[nodiscard]] __device__ DeviceVector& getGrads() noexcept { return grads; }
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

#include "DiffVectorImpl/DiffVectorImpl.cuh"
