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

#include "DiffVectorExpr.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/VectorExpr.cuh"

namespace Physica::Core {
    template<ExpressionType type, unsigned int Order, class T1> class DiffVectorExprUnitary;
    template<ExpressionType type, unsigned int Order, class T1, class T2 = T1> class DiffVectorExprBinary;

    namespace Internal {
        template<ExpressionType Type, class VectorType>
        __global__ void DiffVectorExpr_unitaryOpKernel(
                Physica::PlainStruct<VectorType> result_,
                Physica::PlainStruct<const VectorType> v_) {
            using ScalarType = typename VectorType::ScalarType;
            using DiffRecord = typename VectorType::DiffRecord;
            const VectorType& v = v_.getDerived();
            VectorType& result = result_.getDerived();
            const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
            const size_t length = result.getLength();
            if (index >= length)
                return;
            result.getRecord(index) = DiffRecord{index, Type};
            result.getOperands()[index] = v.calc(index);
            if constexpr (Type == ExpressionType::Relu)
                result.getValue(index) = relu(v.getValue(index));
            else if constexpr (Type == ExpressionType::Exp)
                result.getValue(index) = exp(v.getValue(index));
            else
                static_assert(Type == ExpressionType::Exp, "[Error]: Not implemented");
            result.getGrad(index) = 0;
        }

        template<ExpressionType Type, class VectorType>
        __global__ void DiffVectorExpr_binaryOpKernel(
                Physica::PlainStruct<VectorType> result_,
                Physica::PlainStruct<const VectorType> v_,
                Physica::PlainStruct<const typename VectorType::ScalarType> s_) {
            using ScalarType = typename VectorType::ScalarType;
            using DiffRecord = typename VectorType::DiffRecord;
            const VectorType& v = v_.getDerived();
            const ScalarType& s = s_.getDerived();
            VectorType& result = result_.getDerived();
            const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
            const size_t length = result.getLength();
            if (index >= length)
                return;
            result.getRecord(index) = DiffRecord{index * 2, Type};
            result.getOperands()[index * 2] = v.calc(index);
            result.getOperands()[index * 2 + 1] = s;
            if constexpr (Type == ExpressionType::Add)
                result.getValue(index) = v.getValue(index) + s.getValue();
            else if constexpr (Type == ExpressionType::Sub)
                result.getValue(index) = v.getValue(index) - s.getValue();
            else if constexpr (Type == ExpressionType::Mul)
                result.getValue(index) = v.getValue(index) * s.getValue();
            else if constexpr (Type == ExpressionType::Div)
                result.getValue(index) = v.getValue(index) / s.getValue();
            else
                static_assert(Type == ExpressionType::Add, "[Error]: Not implemented");
            result.getGrad(index) = 0;
        }

        template<ExpressionType Type, class VectorType>
        __global__ void DiffVectorExpr_binaryOpKernel(
                Physica::PlainStruct<VectorType> result_,
                Physica::PlainStruct<const VectorType> v1_,
                Physica::PlainStruct<const VectorType> v2_) {
            using ScalarType = typename VectorType::ScalarType;
            using DiffRecord = typename VectorType::DiffRecord;
            const VectorType& v1 = v1_.getDerived();
            const VectorType& v2 = v2_.getDerived();
            VectorType& result = result_.getDerived();
            const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
            const size_t length = result.getLength();
            if (index >= length)
                return;
            result.getRecord(index) = DiffRecord{index * 2, Type};
            result.getOperands()[index * 2] = v1.calc(index);
            result.getOperands()[index * 2 + 1] = v2.calc(index);
            if constexpr (Type == ExpressionType::Add)
                result.getValue(index) = v1.getValue(index) + v2.getValue(index);
            else if constexpr (Type == ExpressionType::Sub)
                result.getValue(index) = v1.getValue(index) - v2.getValue(index);
            else if constexpr (Type == ExpressionType::Mul)
                result.getValue(index) = v1.getValue(index) * v2.getValue(index);
            else if constexpr (Type == ExpressionType::Div)
                result.getValue(index) = v1.getValue(index) / v2.getValue(index);
            else
                static_assert(Type == ExpressionType::Add, "[Error]: Not implemented");
            result.getGrad(index) = 0;
        }
    }

    template<ExpressionType Type, unsigned int Order, class PlainScalar>
    class DiffVectorExprUnitary<Type, Order, Vector<PlainScalar>> {
        using VectorType = device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse, Order>>;
    public:
        [[nodiscard]] static VectorType calc(const VectorType& v) {
            const size_t length = v.getLength();
            assert(length > 0 && "[Error]: A empty vector does nothing");
            VectorType result(length, Type);
            const size_t numThread = std::min(length, VectorType::MaxThreadPerBlock);
            const size_t numBlock = (length + numThread - 1) / numThread;
            Internal::DiffVectorExpr_unitaryOpKernel<Type, VectorType>
                    <<<numBlock, numThread, 0, StreamPool::getStream()>>>(asStruct(result), asStruct(v));
            return result;
        }
    };

    template<ExpressionType Type, unsigned int Order, class PlainScalar>
    class DiffVectorExprBinary<Type, Order, Vector<PlainScalar>, PlainScalar> {
        using VectorType = device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse, Order>>;
        using ScalarType = typename VectorType::ScalarType;
    public:
        [[nodiscard]] static VectorType calc(const VectorType& v, const ScalarType& s) {
            const size_t length = v.getLength();
            assert(length > 0 && "[Error]: A empty vector does nothing");
            VectorType result(length, Type);
            const size_t numThread = std::min(length, VectorType::MaxThreadPerBlock);
            const size_t numBlock = (length + numThread - 1) / numThread;
            Internal::DiffVectorExpr_binaryOpKernel<Type, VectorType>
                    <<<numBlock, numThread, 0, StreamPool::getStream()>>>(asStruct(result), asStruct(v), asStruct(s));
            return result;
        }
    };

    template<ExpressionType Type, unsigned int Order, class PlainScalar>
    class DiffVectorExprBinary<Type, Order, Vector<PlainScalar>> {
        using VectorType = device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse, Order>>;
    public:
        [[nodiscard]] static VectorType calc(const VectorType& v1, const VectorType& v2) {
            const size_t length = v1.getLength();
            assert(length > 0 && "[Error]: A empty vector does nothing");
            assert(length == v2.getLength() && "[Error]: Vector lengths do not match");
            VectorType result(length, Type);
            const size_t numThread = std::min(length, VectorType::MaxThreadPerBlock);
            const size_t numBlock = (length + numThread - 1) / numThread;
            Internal::DiffVectorExpr_binaryOpKernel<Type, VectorType>
                    <<<numBlock, numThread, 0, StreamPool::getStream()>>>(asStruct(result), asStruct(v1), asStruct(v2));
            return result;
        }
    };
    //////////////////////////////////////Operators//////////////////////////////////////
    //////////////////////////////////////Minus//////////////////////////////////////
    //////////////////////////////////////Add//////////////////////////////////////
    template<class PlainScalar, unsigned int Order>
    inline auto operator+(
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse, Order>>& v,
            const device_obj<Differentiable<PlainScalar, DiffMode::Reverse, Order>>& s) {
        return DiffVectorExprBinary<ExpressionType::Add, Order, Vector<PlainScalar>, PlainScalar>::calc(v, s);
    }

    template<class PlainScalar, unsigned int Order>
    inline auto operator+(
            const device_obj<Differentiable<PlainScalar, DiffMode::Reverse, Order>>& s,
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse, Order>>& v) {
        return v + s;
    }

    template<class PlainScalar, unsigned int Order>
    inline auto operator+(
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse, Order>>& v1,
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse, Order>>& v2) {
        return DiffVectorExprBinary<ExpressionType::Add, Order, Vector<PlainScalar>, Vector<PlainScalar>>::calc(v1, v2);
    }
    //////////////////////////////////////Sub//////////////////////////////////////
    template<class PlainScalar, unsigned int Order>
    inline auto operator-(
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse, Order>>& v,
            const device_obj<Differentiable<PlainScalar, DiffMode::Reverse, Order>>& s) {
        return DiffVectorExprBinary<ExpressionType::Sub, Order, Vector<PlainScalar>, PlainScalar>::calc(v, s);
    }

    template<class PlainScalar, unsigned int Order>
    inline auto operator-(
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse, Order>>& v1,
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse, Order>>& v2) {
        return DiffVectorExprBinary<ExpressionType::Sub, Order, Vector<PlainScalar>, Vector<PlainScalar>>::calc(v1, v2);
    }
    //////////////////////////////////////Mul//////////////////////////////////////
    template<class PlainScalar, unsigned int Order>
    inline auto operator*(
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse, Order>>& v,
            const device_obj<Differentiable<PlainScalar, DiffMode::Reverse, Order>>& s) {
        return DiffVectorExprBinary<ExpressionType::Mul, Order, Vector<PlainScalar>, PlainScalar>::calc(v, s);
    }

    template<class PlainScalar, unsigned int Order>
    inline auto operator*(
            const device_obj<Differentiable<PlainScalar, DiffMode::Reverse, Order>>& s,
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse, Order>>& v) {
        return v * s;
    }

    template<class PlainScalar, unsigned int Order>
    inline auto hadamard(
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse, Order>>& v1,
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse, Order>>& v2) {
        return DiffVectorExprBinary<ExpressionType::Mul, Order, Vector<PlainScalar>, Vector<PlainScalar>>::calc(v1, v2);
    }
    //////////////////////////////////////Div//////////////////////////////////////
    template<class PlainScalar, unsigned int Order>
    inline auto operator/(
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse, Order>>& v,
            const device_obj<Differentiable<PlainScalar, DiffMode::Reverse, Order>>& s) {
        return DiffVectorExprBinary<ExpressionType::Div, Order, Vector<PlainScalar>, PlainScalar>::calc(v, s);
    }

    template<class PlainScalar, unsigned int Order>
    inline auto devide(
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse, Order>>& v1,
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse, Order>>& v2) {
        return DiffVectorExprBinary<ExpressionType::Div, Order, Vector<PlainScalar>, Vector<PlainScalar>>::calc(v1, v2);
    }
    //////////////////////////////////////Compare//////////////////////////////////////
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class PlainScalar, unsigned int Order>
    auto relu(const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse, Order>>& v) {
        return DiffVectorExprUnitary<ExpressionType::Relu, Order, Vector<PlainScalar>>::calc(v);
    }

    template<class PlainScalar, unsigned int Order>
    auto exp(const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse, Order>>& v) {
        return DiffVectorExprUnitary<ExpressionType::Exp, Order, Vector<PlainScalar>>::calc(v);
    }
}
