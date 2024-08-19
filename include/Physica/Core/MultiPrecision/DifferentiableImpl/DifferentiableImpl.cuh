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
        template<class ScalarType, unsigned int Order>
        __global__ void __launch_bounds__(1, 1) Differentiable_minusKernel(
                Physica::PlainStruct<device_obj<TraceSegment<ScalarType, Order>>> segment_,
                Physica::PlainStruct<const device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>> s_) {
            using SegmentType = device_obj<TraceSegment<ScalarType, Order>>;
            using DiffRecord = typename SegmentType::DiffRecord;
            auto& segment = segment_.getDerived();
            segment.getRecords()[0] = DiffRecord{0, ExpressionType::Minus};
            const auto& s = s_.getDerived();
            segment.getOperands()[0] = s;
            segment.getValues()[0] = -s.getValue();
            segment.getGrads()[0] = 0;
        }

        template<class ScalarType, unsigned int Order, ExpressionType Type>
        __global__ void __launch_bounds__(1, 1) Differentiable_calcKernel(
                Physica::PlainStruct<device_obj<TraceSegment<ScalarType, Order>>> segment_,
                Physica::PlainStruct<const device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>> s1_,
                Physica::PlainStruct<const device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>> s2_) {
            using SegmentType = device_obj<TraceSegment<ScalarType, Order>>;
            using DiffRecord = typename SegmentType::DiffRecord;
            auto& segment = segment_.getDerived();
            segment.getRecords()[0] = DiffRecord{0, Type};
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
                static_assert(Type == ExpressionType::Add, "[Error]: Not implemented");
            segment.getGrads()[0] = 0;
        }

        template<class ScalarType, unsigned int Order>
        __global__ void __launch_bounds__(1, 1) Differentiable_reverseKernel(
                Physica::PlainStruct<device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>> s) {
            s.getDerived().setGrad(1);
        }
    }

    template<class ScalarType, unsigned int Order>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>::device_obj(ScalarType s)
            : device_obj(TracerType::getInstance().pushSegment(s)[0]) {}

    template<class ScalarType, unsigned int Order>
    __host__ __device__ device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>::device_obj(
            const ScalarType* pValue_, const ScalarType* pGrad_)
            : pValue(const_cast<ScalarType*>(pValue_)), pGrad(const_cast<ScalarType*>(pGrad_)) {}

    template<class ScalarType, unsigned int Order>
    __host__ __device__ inline bool device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>::operator==(const This& other) const {
        const bool result = pValue == other.pValue;
        assert(result == (pGrad == other.pGrad) && "[Error]: Bad scalar");
        return result;
    }

    template<class ScalarType, unsigned int Order>
    inline device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>
    device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>::operator-() const {
        auto& segment = TracerType::getInstance().pushSegment(1, ExpressionType::Minus);
        Internal::Differentiable_minusKernel<ScalarType, Order><<<1, 1, 0, CUDAContext::getInstance()>>>(asStruct(segment), asStruct(*this));
        return segment[0];
    }

    template<class ScalarType, unsigned int Order>
    inline void device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>::reverse() {
        Internal::Differentiable_reverseKernel<ScalarType, Order><<<1, 1, 0, CUDAContext::getInstance()>>>(asStruct(*this));
        TracerType::getInstance().reverse_from(*this);
    }

    template<class ScalarType, unsigned int Order>
    inline void device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>::reverse_to(This to) {
        Internal::Differentiable_reverseKernel<ScalarType, Order><<<1, 1, 0, CUDAContext::getInstance()>>>(asStruct(*this));
        TracerType::getInstance().reverse(*this, to);
    }

    template<class ScalarType, unsigned int Order>
    ScalarType device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>::toHost_value() const {
        ScalarType result;
        toHostAsync_value(result);
        CUDAContext::getInstance().wait();
        return result;
    }

    template<class ScalarType, unsigned int Order>
    ScalarType device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>::toHost_grad() const {
        ScalarType result;
        toHostAsync_grad(result);
        CUDAContext::getInstance().wait();
        return result;
    }

    template<class ScalarType, unsigned int Order>
    void device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>::toHostAsync_value(ScalarType& value) const {
        check(cudaMemcpyAsync(&value, pValue, sizeof(ScalarType), cudaMemcpyKind::cudaMemcpyDeviceToHost, CUDAContext::getInstance()));
    }

    template<class ScalarType, unsigned int Order>
    void device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>::toHostAsync_grad(ScalarType& grad) const {
        check(cudaMemcpyAsync(&grad, pGrad, sizeof(ScalarType), cudaMemcpyKind::cudaMemcpyDeviceToHost, CUDAContext::getInstance()));
    }

    template<class ScalarType, unsigned int Order>
    __host__ __device__ inline void device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>::swap(device_obj& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(pValue, obj.pValue);
        std::swap(pGrad, obj.pGrad);
    }
    ////////////////////////////////////////////////////////////
    template<class ScalarType, unsigned int Order>
    [[nodiscard]] inline device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>
    operator+(const device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>& s1,
              const device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>& s2) {
        using TracerType = device_obj<DiffTracer<ScalarType, Order>>;
        auto& segment = TracerType::getInstance().pushSegment(1, ExpressionType::Add);
        Internal::Differentiable_calcKernel<ScalarType, Order, ExpressionType::Add>
                <<<1, 1, 0, CUDAContext::getInstance()>>>(asStruct(segment), asStruct(s1), asStruct(s2));
        return segment[0];
    }

    template<class ScalarType, unsigned int Order>
    [[nodiscard]] inline device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>
    operator-(const device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>& s1,
              const device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>& s2) {
        using TracerType = device_obj<DiffTracer<ScalarType, Order>>;
        auto& segment = TracerType::getInstance().pushSegment(1, ExpressionType::Sub);
        Internal::Differentiable_calcKernel<ScalarType, Order, ExpressionType::Sub>
                <<<1, 1, 0, CUDAContext::getInstance()>>>(asStruct(segment), asStruct(s1), asStruct(s2));
        return segment[0];
    }

    template<class ScalarType, unsigned int Order>
    [[nodiscard]] inline device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>
    operator*(const device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>& s1,
              const device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>& s2) {
        using TracerType = device_obj<DiffTracer<ScalarType, Order>>;
        auto& segment = TracerType::getInstance().pushSegment(1, ExpressionType::Mul);
        Internal::Differentiable_calcKernel<ScalarType, Order, ExpressionType::Mul>
                <<<1, 1, 0, CUDAContext::getInstance()>>>(asStruct(segment), asStruct(s1), asStruct(s2));
        return segment[0];
    }

    template<class ScalarType, unsigned int Order>
    [[nodiscard]] inline device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>
    operator/(const device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>& s1,
              const device_obj<Differentiable<ScalarType, DiffMode::Reverse, Order>>& s2) {
        using TracerType = device_obj<DiffTracer<ScalarType, Order>>;
        auto& segment = TracerType::getInstance().pushSegment(1, ExpressionType::Div);
        Internal::Differentiable_calcKernel<ScalarType, Order, ExpressionType::Div>
                <<<1, 1, 0, CUDAContext::getInstance()>>>(asStruct(segment), asStruct(s1), asStruct(s2));
        return segment[0];
    }
}
