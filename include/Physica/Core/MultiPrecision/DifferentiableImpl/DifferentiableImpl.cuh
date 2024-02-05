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
        template<class ScalarType, ExpressionType Type>
        __global__ void __launch_bound__(1, 1) Differentiable_calcKernel(
                Physica::PlainStruct<device_obj<TraceSegment<ScalarType>>> segment_,
                Physica::PlainStruct<const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>> s1_,
                Physica::PlainStruct<const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>> s2_>) {
            using SegmentType = typename decltype(s)::Derived;
            using DiffRecord = typename SegmentType::DiffRecord;
            auto& segment = segment_.getDerived();
            segment.getRecords()[0] = DiffRecord(0, Type);
            const auto& s1 = s1_.getDerived();
            const auto& s2 = s2_.getDerived();
            segment.getOperands()[0] = s1;
            segment.getOperands()[1] = s2;
            if constexpr (Type == ExpressionType::Add)
                segment.getValues()[0] = s1.getValue() + s2.getValue();
            else if constexpr (Type == ExpressionType::Sub)
                segment.getValues()[0] = s1.getValue() - s2.getValue();
            else if constexpr (Type == ExpressionType::Mul)
                segment.getValues()[0] = s1.getValue() * s2.getValue();
            else if constexpr (Type == ExpressionType::Div)
                segment.getValues()[0] = s1.getValue() / s2.getValue();
            else
                static_assert(false, "[Error]: Not implemented");
            segment.getGrads()[0] = 0;
        }
    }

    template<class ScalarType>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse>>::device_obj(ScalarType s)
            : device_obj(TracerType::getInstance().pushSegment(s)[0]) {}

    template<class ScalarType>
    __host__ __device__ device_obj<Differentiable<ScalarType, DiffMode::Reverse>>::device_obj(ScalarType* pValue_, ScalarType* pGrad_)
            : pValue(pValue_), pGrad(pGrad_) {}

    template<class ScalarType>
    __host__ __device__ inline bool device_obj<Differentiable<ScalarType, DiffMode::Reverse>>::operator==(const This& other) const {
        const bool result = pValue == other.pValue;
        assert(result == (pGrad == other.pGrad) && "[Error]: Bad scalar");
        return result;
    }

    template<class ScalarType>
    void device_obj<Differentiable<ScalarType, DiffMode::Reverse>>::toHostAsync_value(ScalarType& value) const {
        cudaCheck(cudaMemcpyAsync(&value, pValue, sizeof(ScalarType), cudaMemcpyKind::cudaMemcpyDeviceToHost, StreamPool::getStream()));
    }

    template<class ScalarType>
    void device_obj<Differentiable<ScalarType, DiffMode::Reverse>>::toHostAsync_grad(ScalarType& grad) const {
        cudaCheck(cudaMemcpyAsync(&grad, pGrad, sizeof(ScalarType), cudaMemcpyKind::cudaMemcpyDeviceToHost, StreamPool::getStream()));
    }

    template<class ScalarType>
    __host__ __device__ inline void device_obj<Differentiable<ScalarType, DiffMode::Reverse>>::swap(device_obj& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(pValue, obj.pValue);
        std::swap(pGrad, obj.pGrad);
    }
    ////////////////////////////////////////////////////////////
    template<class ScalarType>
    [[nodiscard]] inline device_obj<Differentiable<ScalarType, DiffMode::Reverse>>
    operator+(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s1,
              const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s2) {
        using TracerType = device_obj<DiffTracer<ScalarType>>;
        auto& segment = TracerType::getInstance().pushSegment(1, ExpressionType::Add);
        Internal::Differentiable_calcKernel<ScalarType, ExpressionType::Add>
                <<<1, 1, 0, StreamPool::getStream()>>>(asStruct(segment), asStruct(s1), asStruct(s2));
        return segment[0];
    }

    template<class ScalarType>
    [[nodiscard]] inline device_obj<Differentiable<ScalarType, DiffMode::Reverse>>
    operator-(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s1,
              const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s2) {
        using TracerType = device_obj<DiffTracer<ScalarType>>;
        auto& segment = TracerType::getInstance().pushSegment(1, ExpressionType::Sub);
        Internal::Differentiable_calcKernel<ScalarType, ExpressionType::Sub>
                <<<1, 1, 0, StreamPool::getStream()>>>(asStruct(segment), asStruct(s1), asStruct(s2));
        return segment[0];
    }

    template<class ScalarType>
    [[nodiscard]] inline device_obj<Differentiable<ScalarType, DiffMode::Reverse>>
    operator*(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s1,
              const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s2) {
        using TracerType = device_obj<DiffTracer<ScalarType>>;
        auto& segment = TracerType::getInstance().pushSegment(1, ExpressionType::Mul);
        Internal::Differentiable_calcKernel<ScalarType, ExpressionType::Mul>
                <<<1, 1, 0, StreamPool::getStream()>>>(asStruct(segment), asStruct(s1), asStruct(s2));
        return segment[0];
    }

    template<class ScalarType>
    [[nodiscard]] inline device_obj<Differentiable<ScalarType, DiffMode::Reverse>>
    operator/(const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s1,
              const device_obj<Differentiable<ScalarType, DiffMode::Reverse>>& s2) {
        using TracerType = device_obj<DiffTracer<ScalarType>>;
        auto& segment = TracerType::getInstance().pushSegment(1, ExpressionType::Div);
        Internal::Differentiable_calcKernel<ScalarType, ExpressionType::Div>
                <<<1, 1, 0, StreamPool::getStream()>>>(asStruct(segment), asStruct(s1), asStruct(s2));
        return segment[0];
    }
}
