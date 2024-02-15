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

    template<class PlainScalar>
    device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::device_obj(size_t length, ExpressionType type)
            : traceSeg(asStruct(TracerType::getInstance().pushSegment(length, type))) {}

    template<class PlainScalar>
    device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::device_obj(const PlainVector& values)
            : traceSeg(asStruct(TracerType::getInstance().pushSegment(values))) {}

    template<class PlainScalar>
    template<class RandomGenerator>
    inline void device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::random_uniform(RandomGenerator& gen) {
        *this = random_uniform(getLength(), gen);
    }

    template<class PlainScalar>
    template<class RandomGenerator>
    inline void device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::random_normal(RandomGenerator& gen) {
        *this = random_normal(getLength(), gen);
    }

    template<class PlainScalar>
    template<class Distribution, class RandomGenerator>
    inline void device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::random_any(Distribution& dist, RandomGenerator& gen) {
        *this = random_any(getLength(), dist, gen);
    }

    template<class PlainScalar>
    __host__ __device__ inline typename device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::ScalarType
    device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::calc(size_t index) const {
        assert(index < getLength() && "[Error]: Index out of range");
        return getTraceSegment()[index];
    }

    template<class PlainScalar>
    device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>
    device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::copy() const {
        This result{};
        const auto& newTrace = TracerType::getInstance().pushSegment(getTraceSegment().copy());
        result.traceSeg = asStruct(newTrace);
        return result;
    }

    template<class PlainScalar>
    template<bool ComputeMax>
    __device__ void device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::minmaxKernelImpl(SegmentType& result) const {
        extern __shared__ PlainScalar buffer[];
        const size_t length = getLength();
        const size_t index = threadIdx.x;

        PlainScalar value = calc(index).getValue();
        int valueIndex = index;
        for (int i = index + blockDim.x; i < length; i += blockDim.x) {
            const PlainScalar temp = calc(i).getValue();
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
            const PlainScalar temp = buffer[index1];
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
        result.getRecords()[0] = {0, ExpressionType::Assign};
        result.getOperands()[0] = trace[valueIndex];
        result.getValues()[0] = value;
        result.getGrads()[0] = 0;
    }

    template<class PlainScalar>
    __device__ void device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::sumKernelImpl(SegmentType& result) const {
        extern __shared__ PlainScalar buffer[];
        const size_t length = getLength();
        const size_t index = threadIdx.x;

        PlainScalar threadSum = 0;
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
        result.getRecords()[0] = {0, ExpressionType::Sum};
        result.getOperands()[0] = trace[0];
        result.getOperands()[1] = trace[length - 1];
        result.getValues()[0] = threadSum;
        result.getGrads()[0] = 0;
    }

    template<class PlainScalar>
    __host__ __device__ inline PlainScalar*
    device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::value_ptr(size_t index) const noexcept {
        return calc(index).value_ptr();
    }

    template<class PlainScalar>
    __host__ __device__ inline PlainScalar*
    device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::grad_ptr(size_t index) const noexcept {
        return calc(index).grad_ptr();
    }

    template<class PlainScalar>
    __device__ inline typename device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::DiffRecord&
    device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::getRecord(size_t index) {
        return getTraceSegment().getRecords()[index];
    }

    template<class PlainScalar>
    __device__ inline PlainScalar&
    device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::getValue(size_t index) {
        return getTraceSegment().getValues()[index];
    }

    template<class PlainScalar>
    __device__ inline const PlainScalar&
    device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::getValue(size_t index) const {
        return getTraceSegment().getValues()[index];
    }

    template<class PlainScalar>
    __device__ inline PlainScalar&
    device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::getGrad(size_t index) {
        return getTraceSegment().getGrads()[index];
    }

    template<class PlainScalar>
    __device__ inline const PlainScalar&
    device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::getGrad(size_t index) const {
        return getTraceSegment().getGrads()[index];
    }

    template<class PlainScalar>
    typename device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::ScalarType
    device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::max() const {
        auto& trace = TracerType::getInstance().pushSegment(1, ExpressionType::Assign);
        const size_t length = getLength();
        const size_t numThread = length > MaxThreadPerBlock ? MaxThreadPerBlock : length;
        Internal::DiffVector_minmaxKernel<This, SegmentType, true>
                <<<1, numThread, numThread * sizeof(PlainScalar), StreamPool::getStream()>>>(asStruct(*this), asStruct(trace));
        return ScalarType(trace.getValues().data(), trace.getGrads().data());
    }

    template<class PlainScalar>
    typename device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::ScalarType
    device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::min() const {
        auto& trace = TracerType::getInstance().pushSegment(1, ExpressionType::Assign);
        const size_t length = getLength();
        const size_t numThread = length > MaxThreadPerBlock ? MaxThreadPerBlock : length;
        Internal::DiffVector_minmaxKernel<This, SegmentType, false>
                <<<1, numThread, numThread * sizeof(PlainScalar), StreamPool::getStream()>>>(asStruct(*this), asStruct(trace));
        return ScalarType(trace.getValues().data(), trace.getGrads().data());
    }

    template<class PlainScalar>
    typename device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::ScalarType
    device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::sum() const {
        auto& trace = TracerType::getInstance().pushSegment(1, ExpressionType::Sum);
        const size_t length = getLength();
        const size_t numThread = length > MaxThreadPerBlock ? MaxThreadPerBlock : length;
        Internal::DiffVector_sumKernel<This, SegmentType>
                <<<1, numThread, numThread * sizeof(PlainScalar), StreamPool::getStream()>>>(asStruct(*this), asStruct(trace));
        return ScalarType(trace.getValues().data(), trace.getGrads().data());
    }

    template<class PlainScalar>
    template<class RandomGenerator>
    inline device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>
    device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::random_uniform(size_t len, RandomGenerator& gen) {
        return This(PlainVector::random_uniform(len, gen));
    }

    template<class PlainScalar>
    template<class RandomGenerator>
    inline device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>
    device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::random_normal(size_t len, RandomGenerator& gen) {
        return This(PlainVector::random_normal(len, gen));
    }

    template<class PlainScalar>
    template<class Distribution, class RandomGenerator>
    inline device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>
    device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::random_any(size_t len, Distribution& dist, RandomGenerator& gen) {
        return This(PlainVector::random_any(len, dist, gen));
    }
}
