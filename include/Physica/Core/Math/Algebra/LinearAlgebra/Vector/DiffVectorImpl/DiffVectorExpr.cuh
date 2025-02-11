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

#include "DiffVectorExpr.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/VectorExpr.cuh"

namespace Physica {
    template<ExprType type, int Order, Vector> class DiffVectorExprUnitary;
    template<ExprType type, int Order, Vector T1, class T2 = T1> class DiffVectorExprBinary;

    namespace Internal {
        template<ExprType Type, Vector T>
        __global__ void DiffVectorExpr_unitaryOpKernel(
                Physica::PlainStruct<T> result_,
                Physica::PlainStruct<const T> v_) {
            using ScalarType = T::ScalarType;
            using DiffRecord = T::DiffRecord;
            const T& v = v_.getDerived();
            T& result = result_.getDerived();
            const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
            const size_t length = result.getLength();
            if (index >= length)
                return;
            result.getRecord(index) = DiffRecord{index, Type};
            result.getOperands()[index] = v.calc(index);
            if constexpr (Type == ExprType::Relu)
                result.value(index) = relu(v.value(index));
            else if constexpr (Type == ExprType::Exp)
                result.value(index) = exp(v.value(index));
            else
                static_assert(Type == ExprType::Exp, "[Error]: Not implemented");
            result.grad(index) = 0;
        }

        template<ExprType Type, Vector T>
        __global__ void DiffVectorExpr_binaryOpKernel(
                Physica::PlainStruct<T> result_,
                Physica::PlainStruct<const T> v_,
                Physica::PlainStruct<const T::ScalarType> s_) {
            using ScalarType = T::ScalarType;
            using DiffRecord = T::DiffRecord;
            const T& v = v_.getDerived();
            const ScalarType& s = s_.getDerived();
            T& result = result_.getDerived();
            const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
            const size_t length = result.getLength();
            if (index >= length)
                return;
            result.getRecord(index) = DiffRecord{index * 2, Type};
            result.getOperands()[index * 2] = v.calc(index);
            result.getOperands()[index * 2 + 1] = s;
            if constexpr (Type == ExprType::Add)
                result.value(index) = v.value(index) + s.value();
            else if constexpr (Type == ExprType::Sub)
                result.value(index) = v.value(index) - s.value();
            else if constexpr (Type == ExprType::Mul)
                result.value(index) = v.value(index) * s.value();
            else if constexpr (Type == ExprType::Div)
                result.value(index) = v.value(index) / s.value();
            else
                static_assert(Type == ExprType::Add, "[Error]: Not implemented");
            result.grad(index) = 0;
        }

        template<ExprType Type, Vector T>
        __global__ void DiffVectorExpr_binaryOpKernel(
                Physica::PlainStruct<T> result_,
                Physica::PlainStruct<const T> v1_,
                Physica::PlainStruct<const T> v2_) {
            using ScalarType = T::ScalarType;
            using DiffRecord = T::DiffRecord;
            const T& v1 = v1_.getDerived();
            const T& v2 = v2_.getDerived();
            T& result = result_.getDerived();
            const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
            const size_t length = result.getLength();
            if (index >= length)
                return;
            result.getRecord(index) = DiffRecord{index * 2, Type};
            result.getOperands()[index * 2] = v1.calc(index);
            result.getOperands()[index * 2 + 1] = v2.calc(index);
            if constexpr (Type == ExprType::Add)
                result.value(index) = v1.value(index) + v2.value(index);
            else if constexpr (Type == ExprType::Sub)
                result.value(index) = v1.value(index) - v2.value(index);
            else if constexpr (Type == ExprType::Mul)
                result.value(index) = v1.value(index) * v2.value(index);
            else if constexpr (Type == ExprType::Div)
                result.value(index) = v1.value(index) / v2.value(index);
            else
                static_assert(Type == ExprType::Add, "[Error]: Not implemented");
            result.grad(index) = 0;
        }
    }

    template<ExprType Type, int Order, Scalar T>
    class DiffVectorExprUnitary<Type, Order, VectorND<T>> {
        using VectorType = device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>;
    public:
        [[nodiscard]] static VectorType calc(const VectorType& v) {
            const size_t length = v.getLength();
            assert(length > 0 && "[Error]: A empty vector does nothing");
            VectorType result(length, Type);
            const size_t numThread = std::min(length, VectorType::MaxThreadPerBlock);
            const size_t numBlock = (length + numThread - 1) / numThread;
            Internal::DiffVectorExpr_unitaryOpKernel<Type, VectorType>
                    <<<numBlock, numThread, 0, CUDAContext::getInstance()>>>(asStruct(result), asStruct(v));
            check(cudaGetLastError());
            return result;
        }
    };

    template<ExprType Type, int Order, Scalar T>
    class DiffVectorExprBinary<Type, Order, VectorND<T>, T> {
        using VectorType = device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>;
        using ScalarType = VectorType::ScalarType;
    public:
        [[nodiscard]] static VectorType calc(const VectorType& v, const ScalarType& s) {
            const size_t length = v.getLength();
            assert(length > 0 && "[Error]: A empty vector does nothing");
            VectorType result(length, Type);
            const size_t numThread = std::min(length, VectorType::MaxThreadPerBlock);
            const size_t numBlock = (length + numThread - 1) / numThread;
            Internal::DiffVectorExpr_binaryOpKernel<Type, VectorType>
                    <<<numBlock, numThread, 0, CUDAContext::getInstance()>>>(asStruct(result), asStruct(v), asStruct(s));
            check(cudaGetLastError());
            return result;
        }
    };

    template<ExprType Type, int Order, Scalar T>
    class DiffVectorExprBinary<Type, Order, VectorND<T>> {
        using VectorType = device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>;
    public:
        [[nodiscard]] static VectorType calc(const VectorType& v1, const VectorType& v2) {
            const size_t length = v1.getLength();
            assert(length > 0 && "[Error]: A empty vector does nothing");
            assert(length == v2.getLength() && "[Error]: Vector lengths do not match");
            VectorType result(length, Type);
            const size_t numThread = std::min(length, VectorType::MaxThreadPerBlock);
            const size_t numBlock = (length + numThread - 1) / numThread;
            Internal::DiffVectorExpr_binaryOpKernel<Type, VectorType>
                    <<<numBlock, numThread, 0, CUDAContext::getInstance()>>>(asStruct(result), asStruct(v1), asStruct(v2));
            check(cudaGetLastError());
            return result;
        }
    };
    //////////////////////////////////////Operators//////////////////////////////////////
    //////////////////////////////////////Minus//////////////////////////////////////
    //////////////////////////////////////Add//////////////////////////////////////
    template<Scalar T, int Order>
    inline auto operator+(
            const device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>& v,
            const device_obj<Diff<T, DiffMode::Reverse, Order>>& s) {
        return DiffVectorExprBinary<ExprType::Add, Order, VectorND<T>, T>::calc(v, s);
    }

    template<Scalar T, int Order>
    inline auto operator+(
            const device_obj<Diff<T, DiffMode::Reverse, Order>>& s,
            const device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>& v) {
        return v + s;
    }

    template<Scalar T, int Order>
    inline auto operator+(
            const device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>& v1,
            const device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>& v2) {
        return DiffVectorExprBinary<ExprType::Add, Order, VectorND<T>, VectorND<T>>::calc(v1, v2);
    }
    //////////////////////////////////////Sub//////////////////////////////////////
    template<Scalar T, int Order>
    inline auto operator-(
            const device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>& v,
            const device_obj<Diff<T, DiffMode::Reverse, Order>>& s) {
        return DiffVectorExprBinary<ExprType::Sub, Order, VectorND<T>, T>::calc(v, s);
    }

    template<Scalar T, int Order>
    inline auto operator-(
            const device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>& v1,
            const device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>& v2) {
        return DiffVectorExprBinary<ExprType::Sub, Order, VectorND<T>, VectorND<T>>::calc(v1, v2);
    }
    //////////////////////////////////////Mul//////////////////////////////////////
    template<Scalar T, int Order>
    inline auto operator*(
            const device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>& v,
            const device_obj<Diff<T, DiffMode::Reverse, Order>>& s) {
        return DiffVectorExprBinary<ExprType::Mul, Order, VectorND<T>, T>::calc(v, s);
    }

    template<Scalar T, int Order>
    inline auto operator*(
            const device_obj<Diff<T, DiffMode::Reverse, Order>>& s,
            const device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>& v) {
        return v * s;
    }

    template<Scalar T, int Order>
    inline auto hadamard(
            const device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>& v1,
            const device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>& v2) {
        return DiffVectorExprBinary<ExprType::Mul, Order, VectorND<T>, VectorND<T>>::calc(v1, v2);
    }
    //////////////////////////////////////Div//////////////////////////////////////
    template<Scalar T, int Order>
    inline auto operator/(
            const device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>& v,
            const device_obj<Diff<T, DiffMode::Reverse, Order>>& s) {
        return DiffVectorExprBinary<ExprType::Div, Order, VectorND<T>, T>::calc(v, s);
    }

    template<Scalar T, int Order>
    inline auto devide(
            const device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>& v1,
            const device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>& v2) {
        return DiffVectorExprBinary<ExprType::Div, Order, VectorND<T>, VectorND<T>>::calc(v1, v2);
    }
    //////////////////////////////////////Compare//////////////////////////////////////
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<Scalar T, int Order>
    auto relu(const device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>& v) {
        return DiffVectorExprUnitary<ExprType::Relu, Order, VectorND<T>>::calc(v);
    }

    template<Scalar T, int Order>
    auto exp(const device_obj<Diff<VectorND<T>, DiffMode::Reverse, Order>>& v) {
        return DiffVectorExprUnitary<ExprType::Exp, Order, VectorND<T>>::calc(v);
    }
}
