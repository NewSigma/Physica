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

#include <Physica/Core/Parallel/CUDAContext.cuh>
#include <Physica/Core/Parallel/Future/StreamFuture.cuh>
#include <Physica/Core/MultiPrecision/Differentiable.cuh>
#include "TraceSegment.h"

namespace Physica::Core {
    template<class ScalarType, unsigned int Order>
    class device_obj<TraceSegment<ScalarType, Order>> {
        static_assert(!ScalarType::isDifferentiable, "[Error]: Differentiable<> pack is not necessary");
        static_assert(Order > 0, "[Error]: 0 order is not differentiable");
        using host_obj = TraceSegment<ScalarType, Order>;
        using This = device_obj<host_obj>;
        using HostDiffScalar = typename host_obj::DiffScalar;
        using HostValueVector = typename host_obj::ValueVector;
        using HostGradVector = typename host_obj::GradVector;
    public:
        using DiffScalar = device_obj<HostDiffScalar>;
        using DiffRecord = typename host_obj::DiffRecord;
        using RecordArray = Utils::device_obj<Utils::Array<DiffRecord>>;
        using OperandArray = Utils::device_obj<Utils::Array<HostDiffScalar>>;
        using ValueVector = device_obj<HostValueVector>;
        using GradVector = device_obj<HostGradVector>;
        constexpr static size_t DefaultSize = host_obj::DefaultSize;
        constexpr static size_t MaxThreadPerBlock = ValueVector::MaxThreadPerBlock;
    private:
        RecordArray records;
        OperandArray operands;
        ValueVector values;
        GradVector grads;
    public:
        device_obj() = default;
        device_obj(size_t size, ExpressionType type);
        explicit device_obj(ScalarType value);
        explicit device_obj(const HostValueVector& values_);
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] __host__ __device__ inline DiffScalar operator[](size_t index);
        [[nodiscard]] __host__ __device__ inline const DiffScalar operator[](size_t index) const;
        /* Operations */
        void reverse();
        void zero_grad();
        [[nodiscard]] This copy() const;
        void swap(device_obj& __restrict obj) noexcept;

        __device__ void reverseKernelImpl();
        __device__ void copyKernelImpl(This& target) const;
        /* Getters */
        [[nodiscard]] __device__ RecordArray& getRecords() noexcept { return records; }
        [[nodiscard]] __device__ OperandArray& getOperands() noexcept { return operands; }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return records.getLength(); }
        [[nodiscard]] __host__ __device__ size_t getCapacity() const noexcept { return records.getCapacity(); }
        [[nodiscard]] __host__ __device__ bool empty() const noexcept { return records.empty(); }
        [[nodiscard]] __host__ __device__ bool full() const noexcept { return records.full(); }
        [[nodiscard]] size_t getNumOperands() const noexcept { return operands.getLength(); }
        [[nodiscard]] __host__ __device__ bool isFound(DiffScalar s) const noexcept { return find(s) < getLength(); }
        [[nodiscard]] __host__ __device__ size_t find(DiffScalar s) const noexcept;
        [[nodiscard]] __host__ __device__ ValueVector& getValues() noexcept { return values; }
        [[nodiscard]] __host__ __device__ const ValueVector& getValues() const noexcept { return values; }
        [[nodiscard]] __host__ __device__ GradVector& getGrads() noexcept { return grads; }
        [[nodiscard]] __host__ __device__ const GradVector& getGrads() const noexcept { return grads; }
        /* Static members */
        [[nodiscard]] constexpr static unsigned int numOperand(ExpressionType type) { return host_obj::numOperand(type); }
    private:
        /* Operations */
        inline void init(ScalarType value);
        void init(const HostValueVector& values_);
        /* Friends */
        friend class device_obj<DiffTracer<ScalarType, Order>>;
    };
}

#include "TraceSegmentImpl.cuh"
