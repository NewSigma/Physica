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

#include "DiffVectorExpression.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/VectorExpression.cuh"

namespace Physica::Core {
    template<ExpressionType type, class T1, class T2 = T1> class DiffVectorExpression;

    namespace Internal {
        template<ExpressionType Type, class VectorType>
        __global__ void DiffVectorExpression_calcKernel(
                Physica::PlainStruct<VectorType> result_,
                Physica::PlainStruct<const VectorType> v_,
                Physica::PlainStruct<const typename VectorType::DiffScalar> s_) {
            using DiffScalar = typename VectorType::DiffScalar;
            using DiffRecord = typename VectorType::DiffRecord;
            const VectorType& v = v_.getDerived();
            const DiffScalar& s = s_.getDerived();
            VectorType& result = result_.getDerived();
            const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
            const size_t length = result.getLength();
            if (index >= length)
                return;
            result.getRecords()[index] = DiffRecord(index * 2, Type);
            result.getOperands()[index] = DiffScalar(v.value_ptr(index * 2), v.grad_ptr(index * 2));
            result.getOperands()[index + 1] = DiffScalar(s.value_ptr(), s.grad_ptr());
            if constexpr (Type == ExpressionType::Add)
                result.getValues()[index] = v[index].getValue() + s.getValue();
            else if constexpr (Type == ExpressionType::Sub)
                result.getValues()[index] = v[index].getValue() - s.getValue();
            else if constexpr (Type == ExpressionType::Mul)
                result.getValues()[index] = v[index].getValue() * s.getValue();
            else if constexpr (Type == ExpressionType::Div)
                result.getValues()[index] = v[index].getValue() / s.getValue();
            else
                static_assert(false, "[Error]: Not implemented");
            result.getGrads()[index] = 0;
        }

        template<ExpressionType Type, class VectorType>
        __global__ void DiffVectorExpression_calcKernel(
                Physica::PlainStruct<VectorType> result_,
                Physica::PlainStruct<const VectorType> v1_,
                Physica::PlainStruct<const VectorType> v2_) {
            using DiffScalar = typename VectorType::DiffScalar;
            using DiffRecord = typename VectorType::DiffRecord;
            const VectorType& v1 = v1_.getDerived();
            const VectorType& v2 = v2_.getDerived();
            VectorType& result = result_.getDerived();
            const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
            const size_t length = result.getLength();
            if (index >= length)
                return;
            result.getRecords()[index] = DiffRecord(index * 2, Type);
            result.getOperands()[index] = DiffScalar(v1.value_ptr(index * 2), v1.grad_ptr(index * 2));
            result.getOperands()[index + 1] = DiffScalar(v2.value_ptr(index * 2 + 1), v2.grad_ptr(index * 2 + 1));
            if constexpr (Type == ExpressionType::Add)
                result.getValues()[index] = v1[index].getValue() + v2[index].getValue();
            else if constexpr (Type == ExpressionType::Sub)
                result.getValues()[index] = v1[index].getValue() - v2[index].getValue();
            else if constexpr (Type == ExpressionType::Mul)
                result.getValues()[index] = v1[index].getValue() * v2[index].getValue();
            else if constexpr (Type == ExpressionType::Div)
                result.getValues()[index] = v1[index].getValue() / v2[index].getValue();
            else
                static_assert(false, "[Error]: Not implemented");
            result.getGrads()[index] = 0;
        }
    }

    template<ExpressionType Type, class PlainScalar>
    class DiffVectorExpression<Type, Vector<PlainScalar>, PlainScalar> {
        using VectorType = device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>;
        using DiffScalar = typename VectorType::DiffScalar;
    public:
        [[nodiscard]] static VectorType calc(const VectorType& v, const DiffScalar& s) {
            const size_t length = v.getLength();
            VectorType result(length, Type);
            const size_t numThread = length > VectorType::MaxThreadPerBlock ? VectorType::MaxThreadPerBlock : length;
            const size_t numBlock = (length + numThread - 1) / numThread;
            Internal::DiffVectorExpression_calcKernel<Type, VectorType>
                    <<<numBlock, numThread, 0, StreamPool::getStream()>>>(asStruct(result), asStruct(v), asStruct(s));
            return result;
        }
    };

    template<ExpressionType Type, class PlainScalar>
    class DiffVectorExpression<Type, Vector<PlainScalar>> {
        using VectorType = device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>;
    public:
        [[nodiscard]] static VectorType calc(const VectorType& v1, const VectorType& v2) {
            const size_t length = v.getLength();
            assert(length == v2.getLength() && "[Error]: Vector lengths do not match");
            VectorType result(length, Type);
            const size_t numThread = length > VectorType::MaxThreadPerBlock ? VectorType::MaxThreadPerBlock : length;
            const size_t numBlock = (length + numThread - 1) / numThread;
            Internal::DiffVectorExpression_calcKernel<Type, VectorType>
                    <<<numBlock, numThread, 0, StreamPool::getStream()>>>(asStruct(result), asStruct(v1), asStruct(v2));
            return result;
        }
    };
    //////////////////////////////////////Operators//////////////////////////////////////
    //////////////////////////////////////Minus//////////////////////////////////////
    //////////////////////////////////////Add//////////////////////////////////////
    template<class PlainScalar>
    inline auto operator+(
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>& v,
            const device_obj<Differentiable<PlainScalar, DiffMode::Reverse>>& s) {
        return DiffVectorExpression<ExpressionType::Add, PlainScalar, ScalarBase<PlainScalar>>::calc(v, s);
    }

    template<class PlainScalar>
    inline auto operator+(
            const device_obj<Differentiable<PlainScalar, DiffMode::Reverse>>& s,
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>& v) {
        return v + s;
    }

    template<class PlainScalar>
    inline auto operator+(
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>& v1,
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>& v2) {
        return DiffVectorExpression<ExpressionType::Add, Vector<PlainScalar>, Vector<PlainScalar>>::calc(v1, v2);
    }
    //////////////////////////////////////Sub//////////////////////////////////////
    template<class PlainScalar>
    inline auto operator-(
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>& v,
            const device_obj<Differentiable<PlainScalar, DiffMode::Reverse>>& s) {
        return DiffVectorExpression<ExpressionType::Sub, Vector<PlainScalar>, ScalarBase<PlainScalar>>::calc(v, s);
    }

    template<class PlainScalar>
    inline auto operator-(
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>& v1,
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>& v2) {
        return DiffVectorExpression<ExpressionType::Sub, Vector<PlainScalar>, Vector<PlainScalar>>::calc(v1, v2);
    }
    //////////////////////////////////////Mul//////////////////////////////////////
    template<class PlainScalar>
    inline auto operator*(
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>& v,
            const device_obj<Differentiable<PlainScalar, DiffMode::Reverse>>& s) {
        return DiffVectorExpression<ExpressionType::Mul, Vector<PlainScalar>, ScalarBase<PlainScalar>>::calc(v, s);
    }

    template<class PlainScalar>
    inline auto operator*(
            const device_obj<Differentiable<PlainScalar, DiffMode::Reverse>>& s,
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>& v) {
        return v * s;
    }

    template<class PlainScalar>
    inline auto hadamard(
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>& v1,
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>& v2) {
        return DiffVectorExpression<ExpressionType::Mul, Vector<PlainScalar>, Vector<PlainScalar>>::calc(v1, v2);
    }
    //////////////////////////////////////Div//////////////////////////////////////
    template<class PlainScalar>
    inline auto operator/(
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>& v,
            const device_obj<Differentiable<PlainScalar, DiffMode::Reverse>>& s) {
        return DiffVectorExpression<ExpressionType::Div, Vector<PlainScalar>, ScalarBase<PlainScalar>>::calc(v, s);
    }

    template<class PlainScalar>
    inline auto devide(
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>& v1,
            const device_obj<Differentiable<Vector<PlainScalar>, DiffMode::Reverse>>& v2) {
        return DiffVectorExpression<ExpressionType::Div, Vector<PlainScalar>, Vector<PlainScalar>>::calc(v1, v2);
    }
    //////////////////////////////////////Compare//////////////////////////////////////
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
}
