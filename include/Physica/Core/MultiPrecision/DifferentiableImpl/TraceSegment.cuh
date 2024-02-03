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

#include "TraceSegment.h"
#include "Physica/Core/Parallel/StreamPool.cuh"
#include "Physica/Core/Parallel/Future/StreamFuture.cuh"
#include "Physica/Core/MultiPrecision/Differentiable.cuh"

namespace Physica::Core {
    template<class ScalarType>
    class device_obj<TraceSegment<ScalarType>> {
        static_assert(!ScalarType::isDifferentiable, "[Error]: Differentiable<> pack is not necessary");
        using host_obj = TraceSegment<ScalarType>;
        using This = device_obj<host_obj>;
        using DiffRecord = typename host_obj::DiffRecord;
        using HostDiffScalar = typename host_obj::DiffScalar;
        using DiffScalar = device_obj<HostDiffScalar>;
        using RecordArray = Utils::device_obj<Utils::Array<DiffRecord>>;
        using OperandArray = Utils::device_obj<Utils::Array<HostDiffScalar>>;
        using VectorType = typename host_obj::VectorType;
        using DeviceVector = device_obj<VectorType>;
    public:
        constexpr static size_t DefaultSize = host_obj::DefaultSize;
    private:
        RecordArray records;
        OperandArray operands;
        DeviceVector values;
        DeviceVector grads;
    public:
        device_obj() = default;
        explicit device_obj(size_t size);
        explicit device_obj(ExpressionType type);
        device_obj(const VectorType& values_);
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] __host__ __device__ inline DiffScalar operator[](size_t index);
        [[nodiscard]] __host__ __device__ inline const DiffScalar operator[](size_t index) const;
        /* Operations */
        [[nodiscard]] This copy() const;
        void swap(device_obj& __restrict obj) noexcept;

        __device__ void copyKernelImpl(This& target) const;
        /* Getters */
        [[nodiscard]] __device__ RecordArray& getRecords() noexcept { return records; }
        [[nodiscard]] __device__ OperandArray& getOperands() noexcept { return operands; }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return records.getLength(); }
        [[nodiscard]] __host__ __device__ size_t getCapacity() const noexcept { return records.getCapacity(); }
        [[nodiscard]] __host__ __device__ bool empty() const noexcept { return records.empty(); }
        [[nodiscard]] __host__ __device__ bool full() const noexcept { return records.full(); }
        [[nodiscard]] __host__ __device__ bool isFound(DiffScalar s) const noexcept { return find(s) < getLength(); }
        [[nodiscard]] __host__ __device__ size_t find(DiffScalar s) const noexcept;
        [[nodiscard]] __host__ __device__ DeviceVector& getValues() noexcept { return values; }
        [[nodiscard]] __host__ __device__ const DeviceVector& getValues() const noexcept { return values; }
        [[nodiscard]] __host__ __device__ DeviceVector& getGrads() noexcept { return grads; }
        [[nodiscard]] __host__ __device__ const DeviceVector& getGrads() const noexcept { return grads; }
    private:
        friend class device_obj<DiffTracer<ScalarType>>;
    };
}

#include "TraceSegmentImpl.cuh"
