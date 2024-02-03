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
                const Physica::PlainStruct<This> v,
                Physica::PlainStruct<SegmentType> resultTrace) {
            v.getDerived().template minmaxKernelImpl<ComputeMax>(resultTrace.getDerived());
        }
    }

    template<class PlainScalar>
    device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::device_obj(size_t length)
            : traceSeg(asStruct(TracerType::getInstance().pushSegment(length))) {}

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
    __device__ void device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::minmaxKernelImpl(SegmentType& resultTrace) const {
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

        for (int i = length / 2; i > 1; i /= 2) {
            __syncthreads();
            if (index < i) {
                const int index1 = index + i;
                const PlainScalar temp = buffer[index1];
                const bool flag = ComputeMax ? (temp > value) : (temp < value);
                if (flag) {
                    value = temp;
                    valueIndex = index1;
                    buffer[index] = value;
                }
            }
            else
                return;
        }

        assert(index == 0 && "[Error]: Other thread should have exited");
        const bool isSpecialCase = blockDim.x % 2 != 0;
        if (isSpecialCase) {
            const int index1 = blockDim.x - 1;
            const PlainScalar temp = buffer[index1];
            const bool flag = ComputeMax ? (temp > value) : (temp < value);
            if (flag) {
                value = temp;
                valueIndex = index1;
            }
        }

        const auto& trace = getTraceSegment();
        resultTrace.getRecords()[0] = {0, ExpressionType::Assign};
        resultTrace.getOperands()[0] = {trace.getValues().data_ptr(valueIndex), trace.getGrads().data_ptr(valueIndex)};
        resultTrace.getValues()[0] = value;
        resultTrace.getGrads()[0] = 0;
    }

    template<class PlainScalar>
    typename device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::ScalarType
    device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::max() const noexcept {
        auto& trace = TracerType::getInstance().pushSegment(ExpressionType::Assign);
        const size_t length = getLength();
        const size_t numThread = length > MaxThreadPerBlock ? MaxThreadPerBlock : length;
        Internal::DiffVector_minmaxKernel<This, SegmentType, true>
                <<<1, numThread, length * sizeof(PlainScalar), StreamPool::getStream()>>>(asStruct(*this), asStruct(trace));
        return ScalarType(trace.getValues().data(), trace.getGrads().data());
    }

    template<class PlainScalar>
    typename device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::ScalarType
    device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>::min() const noexcept {
        auto& trace = TracerType::getInstance().pushSegment(ExpressionType::Assign);
        const size_t length = getLength();
        const size_t numThread = length > MaxThreadPerBlock ? MaxThreadPerBlock : length;
        Internal::DiffVector_minmaxKernel<This, SegmentType, false>
                <<<1, numThread, length * sizeof(PlainScalar), StreamPool::getStream()>>>(asStruct(*this), asStruct(trace));
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
