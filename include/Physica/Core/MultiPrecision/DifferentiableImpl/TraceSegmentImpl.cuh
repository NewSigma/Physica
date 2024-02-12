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
        __global__ void __launch_bounds__(1, 1) TraceSegment_setKernel(
                Physica::PlainStruct<Core::device_obj<TraceSegment<ScalarType>>> segment_, ScalarType value) {
            using DiffRecord = typename Core::device_obj<TraceSegment<ScalarType>>::DiffRecord;
            auto& segment = segment_.getDerived();
            segment.getRecords()[0] = DiffRecord{0, ExpressionType::Set};
            segment.getValues()[0] = value;
            segment.getGrads()[0] = 0;
        }

        template<class ScalarType>
        __global__ void TraceSegment_copyKernel(
                Physica::PlainStruct<const Core::device_obj<TraceSegment<ScalarType>>> source,
                Physica::PlainStruct<Core::device_obj<TraceSegment<ScalarType>>> target) {
            source.getDerived().copyKernelImpl(target.getDerived());
        }
    }

    template<class ScalarType>
    device_obj<TraceSegment<ScalarType>>::device_obj(size_t size, ExpressionType type)
            : records(size)
            , operands(host_obj::numOperand(type) * size)
            , values(size)
            , grads(size) {}

    template<class ScalarType>
    device_obj<TraceSegment<ScalarType>>::device_obj(ScalarType value)
            : device_obj(1, ExpressionType::Set) {
        Internal::TraceSegment_setKernel<<<1, 1, 0, StreamPool::getStream()>>>(asStruct(*this), value);
    }
    /**
     * Optimize: Allocate \param values_ using pinned memory shall improve performance
     */
    template<class ScalarType>
    device_obj<TraceSegment<ScalarType>>::device_obj(const VectorType& values_)
            : device_obj(values_.getLength(), ExpressionType::Set) {
        values_.toDeviceAsync(values);
        auto future = StreamFuture::makeFuture();
        cudaMemsetAsync((void*)records.data(), 0, getLength() * sizeof(DiffRecord), StreamPool::getStream());
        zero_grad();
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
    void device_obj<TraceSegment<ScalarType>>::zero_grad() {
        cudaMemsetAsync((void*)grads.data(), 0, getLength() * sizeof(ScalarType), StreamPool::getStream());
    }

    template<class ScalarType>
    device_obj<TraceSegment<ScalarType>> device_obj<TraceSegment<ScalarType>>::copy() const {
        constexpr size_t MaxThreadPerBlock = 256;
        This result(getLength());
        const size_t numThread = getLength() < MaxThreadPerBlock ? getLength() : MaxThreadPerBlock;
        const size_t numBlock = (getLength() + numThread - 1) / numThread;
        Internal::TraceSegment_copyKernel<<<numThread, numBlock, 0, StreamPool::getStream()>>>(asStruct(*this), asStruct(result));
        result.values = values;
        zero_grad();
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
}
