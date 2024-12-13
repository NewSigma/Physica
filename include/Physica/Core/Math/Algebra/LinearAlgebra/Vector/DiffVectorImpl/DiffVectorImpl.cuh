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

namespace Physica::Core {
    namespace Internal {
        template<class This, class SegmentType, bool ComputeMax>
        __global__ void DiffVector_minmaxKernel(
                const Physica::PlainStruct<const This> v,
                Physica::PlainStruct<SegmentType> result) {
            v.getDerived().template minmaxKernelImpl<ComputeMax>(result.getDerived());
        }

        template<class This, class SegmentType>
        __global__ void DiffVector_sumKernel(
                const Physica::PlainStruct<const This> v,
                Physica::PlainStruct<SegmentType> result) {
            v.getDerived().sumKernelImpl(result.getDerived());
        }
    }

    template<Scalar T, int Order>
    device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::device_obj(size_t length, ExprType type)
            : traceSeg(asStruct(TracerType::getInstance().pushSegment(length, type))) {}

    template<Scalar T, int Order>
    device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::device_obj(const PlainVector& values)
            : traceSeg(asStruct(TracerType::getInstance().pushSegment(values))) {}

    template<Scalar T, int Order>
    template<RandomGenerator R>
    inline void device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::random_uniform() {
        *this = random_uniform<R>(getLength());
    }

    template<Scalar T, int Order>
    template<RandomGenerator R>
    inline void device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::random_normal() {
        *this = random_normal<R>(getLength());
    }

    template<Scalar T, int Order>
    template<class Distribution, RandomGenerator R>
    inline void device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::random_any(Distribution& dist) {
        *this = random_any<R>(getLength(), dist);
    }

    template<Scalar T, int Order>
    template<Side Owner>
    __host__ __device__ inline device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::ScalarType
    device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::calc(size_t index) const {
        assert(index < getLength() && "[Error]: Index out of range");
        return getTraceSegment()[index];
    }

    template<Scalar T, int Order>
    template<bool ComputeMax>
    __device__ void device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::minmaxKernelImpl(SegmentType& result) const {
        extern __shared__ T buffer[];
        const size_t length = getLength();
        const size_t index = threadIdx.x;

        T value = calc(index).getValue();
        int valueIndex = index;
        for (int i = index + blockDim.x; i < length; i += blockDim.x) {
            const T temp = calc(i).getValue();
            const bool flag = ComputeMax ? (temp > value) : (temp < value);
            if (flag) {
                value = temp;
                valueIndex = i;
            }
        }
        buffer[index] = value;

        int i0 = blockDim.x;
        int i = (i0 + 1) / 2;
        while (index + i < i0) {
            __syncthreads();
            const int index1 = index + i;
            const T temp = buffer[index1];
            const bool flag = ComputeMax ? (temp > value) : (temp < value);
            if (flag) {
                value = temp;
                valueIndex = index1;
                buffer[index] = value;
            }
            i0 = i;
            i = (i0 + 1) / 2;
        }

        if (index != 0)
            return;
        const auto& trace = getTraceSegment();
        result.getRecords()[0] = {0, ExprType::Assign};
        result.getOperands()[0] = trace[valueIndex];
        result.getValues()[0] = value;
        result.getGrads()[0] = 0;
    }

    template<Scalar T, int Order>
    __device__ void device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::sumKernelImpl(SegmentType& result) const {
        extern __shared__ T buffer[];
        const size_t length = getLength();
        const size_t index = threadIdx.x;

        T threadSum = 0;
        for (int i = index; i < length; i += blockDim.x)
            threadSum += calc(i).getValue();
        buffer[index] = threadSum;

        int i0 = blockDim.x;
        int i = (i0 + 1) / 2;
        while (index + i < i0) {
            __syncthreads();
            threadSum += buffer[index + i];
            buffer[index] = threadSum;
            i0 = i;
            i = (i0 + 1) / 2;
        }

        if (index != 0)
            return;
        const auto& trace = getTraceSegment();
        result.getRecords()[0] = {0, ExprType::Sum};
        result.getOperands()[0] = trace[0];
        result.getOperands()[1] = trace[length - 1];
        result.getValues()[0] = threadSum;
        result.getGrads()[0] = 0;
    }

    template<Scalar T, int Order>
    __host__ __device__ inline T*
    device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::value_ptr(size_t index) const noexcept {
        return calc(index).value_ptr();
    }

    template<Scalar T, int Order>
    __host__ __device__ inline T*
    device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::grad_ptr(size_t index) const noexcept {
        return calc(index).grad_ptr();
    }

    template<Scalar T, int Order>
    __device__ inline device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::DiffRecord&
    device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::getRecord(size_t index) {
        return getTraceSegment().getRecords()[index];
    }

    template<Scalar T, int Order>
    __device__ inline T&
    device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::getValue(size_t index) {
        return getTraceSegment().getValues()[index];
    }

    template<Scalar T, int Order>
    __device__ inline const T&
    device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::getValue(size_t index) const {
        return getTraceSegment().getValues()[index];
    }

    template<Scalar T, int Order>
    __device__ inline T&
    device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::getGrad(size_t index) {
        return getTraceSegment().getGrads()[index];
    }

    template<Scalar T, int Order>
    __device__ inline const T&
    device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::getGrad(size_t index) const {
        return getTraceSegment().getGrads()[index];
    }

    template<Scalar T, int Order>
    device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::ScalarType
    device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::max() const {
        auto& trace = TracerType::getInstance().pushSegment(1, ExprType::Assign);
        const size_t length = getLength();
        const size_t numThread = length > MaxThreadPerBlock ? MaxThreadPerBlock : length;
        Internal::DiffVector_minmaxKernel<This, SegmentType, true>
                <<<1, numThread, numThread * sizeof(T), CUDAContext::getInstance()>>>(asStruct(*this), asStruct(trace));
        check(cudaGetLastError());
        return ScalarType(trace.getValues().data(), trace.getGrads().data());
    }

    template<Scalar T, int Order>
    device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::ScalarType
    device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::min() const {
        auto& trace = TracerType::getInstance().pushSegment(1, ExprType::Assign);
        const size_t length = getLength();
        const size_t numThread = length > MaxThreadPerBlock ? MaxThreadPerBlock : length;
        Internal::DiffVector_minmaxKernel<This, SegmentType, false>
                <<<1, numThread, numThread * sizeof(T), CUDAContext::getInstance()>>>(asStruct(*this), asStruct(trace));
        check(cudaGetLastError());
        return ScalarType(trace.getValues().data(), trace.getGrads().data());
    }

    template<Scalar T, int Order>
    device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::ScalarType
    device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::sum() const {
        auto& trace = TracerType::getInstance().pushSegment(1, ExprType::Sum);
        const size_t length = getLength();
        const size_t numThread = length > MaxThreadPerBlock ? MaxThreadPerBlock : length;
        Internal::DiffVector_sumKernel<This, SegmentType>
                <<<1, numThread, numThread * sizeof(T), CUDAContext::getInstance()>>>(asStruct(*this), asStruct(trace));
        check(cudaGetLastError());
        return ScalarType(trace.getValues().data(), trace.getGrads().data());
    }

    template<Scalar T, int Order>
    template<RandomGenerator R>
    inline device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>
    device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::random_uniform(size_t len) {
        return This(PlainVector::random_uniform<R>(len));
    }

    template<Scalar T, int Order>
    template<RandomGenerator R>
    inline device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>
    device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::random_normal(size_t len) {
        return This(PlainVector::random_normal<R>(len));
    }

    template<Scalar T, int Order>
    template<class Distribution, RandomGenerator R>
    inline device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>
    device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>::random_any(size_t len, Distribution& dist) {
        return This(PlainVector::random_any<R>(len, dist));
    }
}
