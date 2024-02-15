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
        __global__ void __launch_bounds__(256, 1)
        TraceSegment_reverseKernel(Physica::PlainStruct<Core::device_obj<TraceSegment<ScalarType>>> segment) {
            static_assert(Core::device_obj<TraceSegment<ScalarType>>::MaxThreadPerBlock == 256,
                    "[Error]: Keep MaxThreadPerBlock and __launch_bounds__ consistent");
            segment.getDerived().reverseKernelImpl();
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
    void device_obj<TraceSegment<ScalarType>>::reverse() {
        assert(!empty());
        const size_t numThread = getLength() < MaxThreadPerBlock ? getLength() : MaxThreadPerBlock;
        const size_t numBlock = (getLength() + numThread - 1) / numThread;
        Internal::TraceSegment_reverseKernel<<<numBlock, numThread, 0, StreamPool::getStream()>>>(asStruct(*this));
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
        Internal::TraceSegment_copyKernel<<<numBlock, numThread, 0, StreamPool::getStream()>>>(asStruct(*this), asStruct(result));
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
    __device__ void device_obj<TraceSegment<ScalarType>>::reverseKernelImpl() {
        const int index = blockIdx.x * blockDim.x + threadIdx.x;
        if (index >= getLength())
            return;

        const DiffRecord record = records[index];
        if (record.source == ExpressionType::Set)
            return;

        const ScalarType grad = grads[index];
        if (grad.isZero())
            return;
            
        const ScalarType value = values[index];
        const size_t idFirstOperand = record.idFirstOperand;
        const ExpressionType source = record.source;
        DiffScalar operandX = operands[idFirstOperand];
        ScalarType& gradX = operandX.getGrad();
        /* Unitary Operations */ {
            switch (source) {
                case ExpressionType::Assign:
                case ExpressionType::Minus:
                    gradX += grad * ScalarType(source == ExpressionType::Assign ? 1.0 : -1.0);
                    return;
                case ExpressionType::Reciprocal:
                    gradX -= grad * square(value);
                    return;
                case ExpressionType::Sqrt:
                    gradX += grad / value * ScalarType(0.5);
                    return;
                case ExpressionType::Cbrt:
                    gradX += grad / (square(value) * ScalarType(3));
                    return;
                case ExpressionType::Abs:
                    gradX += operandX.getValue().isPositive() ? grad : -grad;
                    return;
                case ExpressionType::Relu:
                    gradX += operandX.getValue().isPositive() ? grad : ScalarType(0);
                    return;
                case ExpressionType::Square:
                    gradX += grad * operandX.getValue() * ScalarType(2);
                    return;
                case ExpressionType::Ln:
                    gradX += grad / operandX.getValue();
                    return;
                case ExpressionType::Exp:
                    gradX += grad * value;
                    return;
                case ExpressionType::Sin:
                    gradX += grad * cos(operandX.getValue());
                    return;
                case ExpressionType::Cos:
                    gradX -= grad * sin(operandX.getValue());
                    return;
                case ExpressionType::ArcCos:
                    gradX -= grad / sqrt(ScalarType(1) - square(operandX.getValue()));
                    return;
                default:;
            }
        }
        /* Binary Operations */ {
            DiffScalar operandY = operands[idFirstOperand + 1];
            ScalarType& gradY = operandY.getGrad();
            switch (source) {
                case ExpressionType::Add:
                case ExpressionType::Sub:
                    gradX += grad;
                    gradY += grad * ScalarType(source == ExpressionType::Add ? 1.0 : -1.0);
                    return;
                case ExpressionType::Mul:
                    gradX += grad * operandY.getValue();
                    gradY += grad * operandX.getValue();
                    return;
                case ExpressionType::Div: {
                    const ScalarType dx = grad * reciprocal(operandY.getValue());
                    gradX += dx;
                    gradY -= dx * value;
                    return;
                }
                case ExpressionType::Sum: {
                    ScalarType* const pGradY = &gradY;
                    for (auto* pGradX = &gradX; pGradX <= pGradY; pGradX += 1)
                        *pGradX += grad;
                    return;
                }
                default:;
            }
        }
        assert(false && "[Error]: Undefined operator for back propagation");
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
