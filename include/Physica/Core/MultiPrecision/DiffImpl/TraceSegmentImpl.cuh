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
        template<Scalar T, int Order>
        __global__ void __launch_bounds__(1, 1) TraceSegment_setKernel(
                Physica::PlainStruct<Core::device_obj<TraceSegment<T, Order>>> segment_, T value) {
            using DiffRecord = Core::device_obj<TraceSegment<T, Order>>::DiffRecord;
            auto& segment = segment_.getDerived();
            segment.getRecords()[0] = DiffRecord{0, ExprType::Set};
            segment.getValues()[0] = value;
            segment.getGrads()[0] = 0;
        }

        template<Scalar T, int Order>
        __global__ void __launch_bounds__(256, 1)
        TraceSegment_reverseKernel(Physica::PlainStruct<Core::device_obj<TraceSegment<T, Order>>> segment) {
            static_assert(Core::device_obj<TraceSegment<T, Order>>::MaxThreadPerBlock == 256,
                    "[Error]: Keep MaxThreadPerBlock and __launch_bounds__ consistent");
            segment.getDerived().reverseKernelImpl();
        }

        template<Scalar T, int Order>
        __global__ void TraceSegment_copyKernel(
                Physica::PlainStruct<const Core::device_obj<TraceSegment<T, Order>>> source,
                Physica::PlainStruct<Core::device_obj<TraceSegment<T, Order>>> target) {
            source.getDerived().copyKernelImpl(target.getDerived());
        }
    }

    template<Scalar T, int Order>
    device_obj<TraceSegment<T, Order>>::device_obj(size_t size, ExprType type)
            : records(size)
            , operands(size * numOperand(type))
            , values(size)
            , grads(size) {}

    template<Scalar T, int Order>
    device_obj<TraceSegment<T, Order>>::device_obj(T value)
            : device_obj(1, ExprType::Set) {
        init(value);
    }
    /**
     * Optimize: Allocate \param values_ using pinned memory shall improve performance
     */
    template<Scalar T, int Order>
    device_obj<TraceSegment<T, Order>>::device_obj(const HostValueVector& values_)
            : device_obj(values_.getLength(), ExprType::Set) {
        init(values_);
    }

    template<Scalar T, int Order>
    __host__ __device__ inline device_obj<TraceSegment<T, Order>>::DiffScalar
    device_obj<TraceSegment<T, Order>>::operator[](size_t index) {
        return DiffScalar(values.data_ptr(index), grads.data_ptr(index));
    }

    template<Scalar T, int Order>
    __host__ __device__ inline const device_obj<TraceSegment<T, Order>>::DiffScalar
    device_obj<TraceSegment<T, Order>>::operator[](size_t index) const {
        return const_cast<This&>(*this).operator[](index);
    }

    template<Scalar T, int Order>
    void device_obj<TraceSegment<T, Order>>::reverse() {
        assert(!empty());
        const size_t numThread = getLength() < MaxThreadPerBlock ? getLength() : MaxThreadPerBlock;
        const size_t numBlock = (getLength() + numThread - 1) / numThread;
        Internal::TraceSegment_reverseKernel<<<numBlock, numThread, 0, CUDAContext::getInstance()>>>(asStruct(*this));
        check(cudaGetLastError());
    }

    template<Scalar T, int Order>
    void device_obj<TraceSegment<T, Order>>::zero_grad() {
        cudaMemsetAsync((void*)grads.data(), 0, getLength() * sizeof(T), CUDAContext::getInstance());
    }

    template<Scalar T, int Order>
    device_obj<TraceSegment<T, Order>> device_obj<TraceSegment<T, Order>>::copy() const {
        constexpr size_t MaxThreadPerBlock = 256;
        This result(getLength());
        const size_t numThread = getLength() < MaxThreadPerBlock ? getLength() : MaxThreadPerBlock;
        const size_t numBlock = (getLength() + numThread - 1) / numThread;
        Internal::TraceSegment_copyKernel<<<numBlock, numThread, 0, CUDAContext::getInstance()>>>(asStruct(*this), asStruct(result));
        check(cudaGetLastError());
        result.values = values;
        zero_grad();
        return result;
    }

    template<Scalar T, int Order>
    TraceSegment<T, Order> device_obj<TraceSegment<T, Order>>::toHost() const {
        host_obj result(getCapacity());
        records.toHost(result.records);
        operands.toHost(result.operands);
        values.toHost(result.values);
        grads.toHost(result.grads);
        return result;
    }

    template<Scalar T, int Order>
    void device_obj<TraceSegment<T, Order>>::swap(device_obj& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        records.swap(obj.records);
        operands.swap(obj.operands);
        values.swap(obj.values);
        grads.swap(obj.grads);
    }

    template<Scalar T, int Order>
    __device__ void device_obj<TraceSegment<T, Order>>::reverseKernelImpl() {
        const int index = blockIdx.x * blockDim.x + threadIdx.x;
        if (index >= getLength())
            return;

        const DiffRecord record = records[index];
        if (record.source == ExprType::Set)
            return;

        const T grad = grads[index];
        if (grad.isZero())
            return;
            
        const T value = values[index];
        const size_t idFirstOperand = record.idFirstOperand;
        const ExprType source = record.source;
        DiffScalar operandX = operands[idFirstOperand];
        T& gradX = operandX.getGrad();
        /* Unitary Operations */ {
            switch (source) {
                case ExprType::Assign:
                case ExprType::Minus:
                    gradX += grad * T(source == ExprType::Assign ? 1.0 : -1.0);
                    return;
                case ExprType::Reciprocal:
                    gradX -= grad * square(value);
                    return;
                case ExprType::Sqrt:
                    gradX += grad / value * T(0.5);
                    return;
                case ExprType::Cbrt:
                    gradX += grad / (square(value) * T(3));
                    return;
                case ExprType::Abs:
                    gradX += operandX.getValue().isPositive() ? grad : -grad;
                    return;
                case ExprType::Relu:
                    gradX += operandX.getValue().isPositive() ? grad : T(0);
                    return;
                case ExprType::Square:
                    gradX += grad * operandX.getValue() * T(2);
                    return;
                case ExprType::Ln:
                    gradX += grad / operandX.getValue();
                    return;
                case ExprType::Exp:
                    gradX += grad * value;
                    return;
                case ExprType::Sin:
                    gradX += grad * cos(operandX.getValue());
                    return;
                case ExprType::Cos:
                    gradX -= grad * sin(operandX.getValue());
                    return;
                case ExprType::ArcCos:
                    gradX -= grad / sqrt(T(1) - square(operandX.getValue()));
                    return;
                default:;
            }
        }
        /* Binary Operations */ {
            DiffScalar operandY = operands[idFirstOperand + 1];
            T& gradY = operandY.getGrad();
            switch (source) {
                case ExprType::Add:
                case ExprType::Sub:
                    gradX += grad;
                    gradY += grad * T(source == ExprType::Add ? 1.0 : -1.0);
                    return;
                case ExprType::Mul:
                    gradX += grad * operandY.getValue();
                    gradY += grad * operandX.getValue();
                    return;
                case ExprType::Div: {
                    const T dx = grad * reciprocal(operandY.getValue());
                    gradX += dx;
                    gradY -= dx * value;
                    return;
                }
                case ExprType::Sum: {
                    T* const pGradY = &gradY;
                    for (auto* pGradX = &gradX; pGradX <= pGradY; pGradX += 1)
                        *pGradX += grad;
                    return;
                }
                default:;
            }
        }
        assert(false && "[Error]: Undefined operator for back propagation");
    }

    template<Scalar T, int Order>
    __device__ void device_obj<TraceSegment<T, Order>>::copyKernelImpl(This& target) const {
        const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
        if (index < getLength()) {
            target.records[index] = DiffRecord{index, ExprType::Assign};
            target.operands[index] = DiffScalar(const_cast<T*>(values.data() + index),
                                                const_cast<T*>(grads.data() + index));
        }
    }

    template<Scalar T, int Order>
    __host__ __device__ size_t device_obj<TraceSegment<T, Order>>::find(DiffScalar s) const noexcept {
        assert(values.getLength() == grads.getLength() && "[Error]: Invalid param");
        const auto* pValue = values.data();
        const auto* pValue1 = s.value_ptr();
        const size_t length = values.getLength();
        if (pValue1 < pValue)
            return length;
        const size_t index = pValue1 - pValue;
        [[maybe_unused]] const bool isValueFound = index < length;
        assert((!isValueFound || (s == this->operator[](index))) && "[Error]: Bad DiffScalar");
        return index;
    }

    template<Scalar T, int Order>
    inline void device_obj<TraceSegment<T, Order>>::init(T value) {
        Internal::TraceSegment_setKernel<<<1, 1, 0, CUDAContext::getInstance()>>>(asStruct(*this), value);
        check(cudaGetLastError());
    }

    template<Scalar T, int Order>
    inline void device_obj<TraceSegment<T, Order>>::init(const HostValueVector& values_) {
        values_.toDeviceAsync(values);
        auto future = StreamFuture::makeFuture();
        cudaMemsetAsync((void*)records.data(), 0, getLength() * sizeof(DiffRecord), CUDAContext::getInstance());
        zero_grad();
        future->wait();
    }
}
