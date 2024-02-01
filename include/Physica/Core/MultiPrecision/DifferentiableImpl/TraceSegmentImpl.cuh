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

namespace Physica::Core {
    namespace Internal {
        template<class ScalarType>
        __global__ void TraceSegment_copyKernel(
                Physica::PlainStruct<Core::device_obj<TraceSegment<ScalarType>>> source,
                Physica::PlainStruct<Core::device_obj<TraceSegment<ScalarType>>> target) {
            source.getDerived().copyKernelImpl(target.getDerived());
        }
    }

    template<class ScalarType>
    device_obj<TraceSegment<ScalarType>>::device_obj(size_t size) {
        records.reserve(size);
        operands.reserve(3 * size); //MulAdd operation is 3-operand
        values.reserve(size);
        grads.reserve(size);
    }
    /**
     * Optimize: Allocate \param values_ using pinned memory shall improve performance
     */
    template<class ScalarType>
    device_obj<TraceSegment<ScalarType>>::device_obj(const VectorType& values_)
            : records(values_.getLength())
            , values(values_.getLength())
            , grads(values_.getLength()) {
        values_.toDeviceAsync(values);
        auto future = StreamFuture::makeFuture();
        cudaMemsetAsync((void*)records.data(), 0, getLength() * sizeof(DiffRecord), StreamPool::getStream());
        cudaMemsetAsync((void*)grads.data(), 0, getLength() * sizeof(DiffScalar), StreamPool::getStream());
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
    device_obj<TraceSegment<ScalarType>> device_obj<TraceSegment<ScalarType>>::copy() const {
        constexpr size_t MaxThreadPerBlock = 256;
        This result(getLength());
        const size_t numThread = getLength() < MaxThreadPerBlock ? getLength() : MaxThreadPerBlock;
        const size_t numBlock = (getLength() + numThread - 1) / numThread;
        Internal::TraceSegment_copyKernel<<<numThread, numBlock, 0, StreamPool::getStream()>>>(asStruct(*this), asStruct(result));
        result.values = values;
        cudaMemsetAsync((void*)grads.data(), 0, getLength() * sizeof(DiffScalar), StreamPool::getStream());
        return result;
    }

    template<class ScalarType>
    void device_obj<TraceSegment<ScalarType>>::swap(device_obj& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        records.swap(obj.records);
        operands.swap(obj.operands);
        values.swap(obj.values);
        grads.swap(obj.grads);
    }

    template<class ScalarType>
    __device__ void device_obj<TraceSegment<ScalarType>>::copyKernelImpl(This& target) const {
        const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
        if (index < getLength()) {
            target.records[index] = DiffRecord{index, ExpressionType::Assign};
            target.operands[index] = DiffScalar(const_cast<ScalarType*>(values.data() + index),
                                                const_cast<ScalarType*>(grads.data() + index));
        }
    }

    template<class ScalarType>
    __host__ __device__ size_t device_obj<TraceSegment<ScalarType>>::find(DiffScalar s) const noexcept {
        return host_obj::findImpl(values, grads, s);
    }

    template<class ScalarType>
    device_obj<TraceSegment<ScalarType>> device_obj<TraceSegment<ScalarType>>::makeSingleNode() {
        This result{};
        result.records.reserve(1);
        result.values.reserve(1);
        result.grads.reserve(1);

        auto future = StreamFuture::makeFuture();
        cudaMemsetAsync((void*)result.records.data(), 0, sizeof(DiffRecord), StreamPool::getStream());
        future->wait();
        return result;
    }
}
