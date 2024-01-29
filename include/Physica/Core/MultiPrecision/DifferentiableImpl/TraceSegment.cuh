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
#include "Physica/Core/MultiPrecision/Differentiable.cuh"

namespace Physica::Core {
    template<class ScalarType>
    class device_obj<TraceSegment<ScalarType>> {
        static_assert(!ScalarType::isDifferentiable, "[Error]: Differentiable<> pack is not necessary");
        using host_obj = TraceSegment<ScalarType>;
        using This = device_obj<This>;
        using DiffRecord = typename host_obj::DiffRecord;
        using DiffScalar = device_obj<typename host_obj::DiffScalar>;
        using VectorType = typename host_obj::VectorType;
        using DeviceVector = device_obj<VectorType>;
        constexpr static size_t DefaultSize = host_obj::DefaultSize;
    private:
        device_obj<Utils::Array<DiffRecord>> records;
        device_obj<Utils::Array<DiffScalar>> operands;
        DeviceVector values;
        DeviceVector grads;
    public:
        explicit device_obj(size_t size = DefaultSize);
        device_obj(const VectorType& values_);
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] __host__ __device__ inline DiffScalar operator[](size_t index);
        [[nodiscard]] __host__ __device__ inline const DiffScalar operator[](size_t index) const;
        /* Operations */
        void swap(device_obj& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return records.getLength(); }
        [[nodiscard]] __host__ __device__ size_t getCapacity() const noexcept { return records.getCapacity(); }
        [[nodiscard]] __host__ __device__ bool empty() const noexcept { return records.empty(); }
        [[nodiscard]] __host__ __device__ bool full() const noexcept { return records.full(); }
    };

    template<class ScalarType>
    device_obj<TraceSegment<ScalarType>>::device_obj(size_t size) {
        assert(size >= DefaultSize && "[Error]: Allocate a small segment maybe bad to performance");
        records.reserve(size);
        operands.reserve(3 * size); //MulAdd operation is 3-operand
        values.reserve(size);
        grads.reserve(size);
    }

    template<class ScalarType>
    device_obj<TraceSegment<ScalarType>>::device_obj(const VectorType& values_)
            : records(values_.getLength())
            , values(values_.getLength())
            , grads(values_.getLength()) {
        values_.toDeviceAsync(values);
        auto future = StreamFuture::makeFuture();
        cudaMemsetAsync(records.data(), 0, getLength() * sizeof(DiffRecord), &StreamPool::getStream());
        cudaMemsetAsync(grads.data_ptr(), 0, getLength() * sizeof(DiffScalar), &StreamPool::getStream());
        future->wait();
    }

    template<class ScalarType>
    __host__ __device__ inline typename device_obj<TraceSegment<ScalarType>>::DiffScalar
    device_obj<TraceSegment<ScalarType>>::operator[](size_t index) {
        return DiffScalar(values.data_ptr(index), grads.data_ptr(index));
    }

    template<class ScalarType>
    __host__ __device__ inline const typename device_obj<TraceSegment<ScalarType>>::DiffScalar
    device_obj<TraceSegment<ScalarType>>::operator[](size_t index) const {
        return DiffScalar(const_cast<ScalarType*>(values.data_ptr(index)), const_cast<ScalarType*>(grads.data_ptr(index)));
    }

    template<class ScalarType>
    void device_obj<TraceSegment<ScalarType>>::swap(device_obj& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        records.swap(obj.records);
        operands.swap(obj.operands);
        values.swap(obj.values);
        grads.swap(obj.grads);
    }
}
