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
        template<class This, class SegmentType>
        __global__ void SGD_stepKernel(Physica::PlainStruct<const This> sgd, Physica::PlainStruct<SegmentType> segment) {
            sgd.getDerived().stepKernelImpl(segment.getDerived());
        }
    }

    template<class ScalarType>
    SGD<device_obj<ScalarType>>::SGD(PlainScalar learnRate_, unsigned int batchSize_) : batchSize(batchSize_) {
        setLearnRate(learnRate_);
    }

    template<class ScalarType>
    inline void SGD<device_obj<ScalarType>>::recordBegin() {
        auto& segment = TracerType::getInstance().pushSegment(1, ExpressionType::Set);
        to = DeviceScalar(segment.getValues().data(), segment.getGrads().data());
    }
    
    template<class ScalarType>
    inline void SGD<device_obj<ScalarType>>::recordEnd() {
        auto& segment = TracerType::getInstance().pushSegment(1, ExpressionType::Set);
        from = DeviceScalar(segment.getValues().data(), segment.getGrads().data());
    }

    template<class ScalarType>
    void SGD<device_obj<ScalarType>>::step() {
        auto& tracer = TracerType::getInstance();
        tracer.forSegmentInRange(from, to, [this](SegmentType& segment) {
            const size_t length = segment.getLength();
            const size_t numThread = length < MaxNumThreadPerBlock ? length : MaxNumThreadPerBlock;
            const size_t numBlock = (length + numThread - 1) / numThread;
            Internal::SGD_stepKernel<This, SegmentType><<<numBlock, numThread, 0, StreamPool::getStream()>>>(asStruct(*this), asStruct(segment));
        });
    }

    template<class ScalarType>
    void SGD<device_obj<ScalarType>>::swap(SGD& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        learnRate.swap(obj.learnRate);
        meanLearnRate.swap(obj.meanLearnRate);
        std::swap(batchSize, obj.batchSize);
        from.swap(obj.from);
        to.swap(obj.to);
    }

    template<class ScalarType>
    __device__ void SGD<device_obj<ScalarType>>::stepKernelImpl(SegmentType& segment) const {
        using DiffScalar = typename SegmentType::DiffScalar;
        const unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;
        if (index < segment.getLength()) {
            DiffScalar s = segment[index];
            s.setValue(s.getValue() - meanLearnRate * s.getGrad());
        }
    }

    template<class ScalarType>
    void SGD<device_obj<ScalarType>>::setLearnRate(PlainScalar lr) {
        learnRate = lr;
        meanLearnRate = lr / PlainScalar(batchSize);
    }
}
