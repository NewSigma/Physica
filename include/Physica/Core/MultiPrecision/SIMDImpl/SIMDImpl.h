/*
 * Copyright 2023 Weibo He.
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

#include "Physica/Core/MultiPrecision/ScalarImpl/ExprType.h"

namespace Physica::Core {
    template<class ScalarType, size_t Size>
    [[nodiscard]] inline ScalarType SIMD<ScalarType, Size>::operator[](int index) const {
        if constexpr (isForward)
            return ScalarType(Base::operator[](index * 2), Base::operator[](index * 2 + 1));
        else
            return ScalarType(Base::operator[](index));
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> SIMD<ScalarType, Size>::operator+(const SIMD& other) const {
        return SIMD(getImpl() + other.getImpl());
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> SIMD<ScalarType, Size>::operator-(const SIMD& other) const {
        return SIMD(getImpl() - other.getImpl());
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> SIMD<ScalarType, Size>::operator*(const SIMD& other) const {
        return SIMD(getImpl() * other.getImpl());
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> SIMD<ScalarType, Size>::operator/(const SIMD& other) const {
        return SIMD(getImpl() / other.getImpl());
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> SIMD<ScalarType, Size>::operator-() const {
        return SIMD(-getImpl());
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline typename SIMD<ScalarType, Size>::BoolSIMDType
    SIMD<ScalarType, Size>::operator>(const SIMD& other) const {
        return BoolSIMDType(getImpl() > other.getImpl());
    }
    
    template<class ScalarType, size_t Size>
    [[nodiscard]] inline typename SIMD<ScalarType, Size>::BoolSIMDType
    SIMD<ScalarType, Size>::operator<(const SIMD& other) const {
        return BoolSIMDType(getImpl() < other.getImpl());
    }

    template<class ScalarType, size_t Size>
    inline void SIMD<ScalarType, Size>::load(const ScalarType* p) {
        Base::load(reinterpret_cast<const TrivialType*>(p));
    }

    template<class ScalarType, size_t Size>
    inline void SIMD<ScalarType, Size>::load_partial(int n, const ScalarType* p) {
        Base::load_partial(n, reinterpret_cast<const TrivialType*>(p));
    }

    template<class ScalarType, size_t Size>
    inline void SIMD<ScalarType, Size>::store(ScalarType* p) const {
        Base::store(reinterpret_cast<TrivialType*>(p));
    }

    template<class ScalarType, size_t Size>
    inline void SIMD<ScalarType, Size>::store_partial(int n, ScalarType* p) const {
        Base::store_partial(n, reinterpret_cast<TrivialType*>(p));
    }

    template<class ScalarType, size_t Size>
    inline void SIMD<ScalarType, Size>::insert(int index, const ScalarType& value) {
        if constexpr (isForward) {
            Base::insert(index * 2, value.getValue().getTrivial());
            Base::insert(index * 2 + 1, value.getGrad().getTrivial());
        }
        else
            Base::insert(index, value.getTrivial());
    }

    template<class ScalarType, size_t Size>
    inline ScalarType SIMD<ScalarType, Size>::horizontal_add() const {
        return Physica::horizontal_add(getImpl());
    }

    template<class ScalarType, size_t Size>
    inline ScalarType SIMD<ScalarType, Size>::horizontal_max() const {
        return Physica::horizontal_max(getImpl());
    }

    template<class ScalarType, size_t Size>
    inline ScalarType SIMD<ScalarType, Size>::horizontal_min() const {
        return Physica::horizontal_min(getImpl());
    }
    //////////////////////////////////////////////////////////////////
    template<class PlainScalar, size_t Size>
    SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size>::SIMD(PlainScalar s) : Base(s) {
        auto& tracer = TracerType::getInstance();
        headNode = tracer.pushOperation(*this, ExprType::Set);
    }

    template<class PlainScalar, size_t Size>
    SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size>::SIMD(ScalarType s) : Base(s.getValue()) {
        auto& tracer = TracerType::getInstance();
        headNode = tracer.pushOperation(*this, ExprType::Assign);
        ScalarType operand[Size];
        for (size_t i = 0; i < Size; ++i)
            operand[i] = s;
        tracer.pushOperand(operand);
    }

    template<class PlainScalar, size_t Size>
    SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size>::SIMD(Base base, ScalarType headNode_)
            : Base(std::move(base)), headNode(headNode_) {}

    template<class PlainScalar, size_t Size>
    inline typename SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size>::ScalarType
    SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size>::operator[](int i) const {
        return ScalarType(value_ptr() + i, grad_ptr() + i);
    }

    template<class PlainScalar, size_t Size>
    inline SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size>
    SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size>::operator+(const SIMD& other) const {
        auto& tracer = TracerType::getInstance();
        const auto temp = Base::operator+(other);
        const ScalarType newHeadNode = tracer.pushOperation(temp, ExprType::Add);
        ScalarType operand[Size * 2];
        for (size_t i = 0; i < Size; ++i) {
            operand[2 * i] = ScalarType(value_ptr() + i, grad_ptr() + i);
            operand[2 * i + 1] = ScalarType(other.value_ptr() + i, other.grad_ptr() + i);
        }
        tracer.pushOperand(operand);
        return {temp, newHeadNode};
    }

    template<class PlainScalar, size_t Size>
    inline SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size>
    SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size>::operator-(const SIMD& other) const {
        auto& tracer = TracerType::getInstance();
        const auto temp = Base::operator-(other);
        const ScalarType newHeadNode = tracer.pushOperation(temp, ExprType::Sub);
        ScalarType operand[Size * 2];
        for (size_t i = 0; i < Size; ++i) {
            operand[2 * i] = ScalarType(value_ptr() + i, grad_ptr() + i);
            operand[2 * i + 1] = ScalarType(other.value_ptr() + i, other.grad_ptr() + i);
        }
        tracer.pushOperand(operand);
        return {temp, newHeadNode};
    }

    template<class PlainScalar, size_t Size>
    inline SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size>
    SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size>::operator*(const SIMD& other) const {
        auto& tracer = TracerType::getInstance();
        const auto temp = Base::operator*(other);
        const ScalarType newHeadNode = tracer.pushOperation(temp, ExprType::Mul);
        ScalarType operand[Size * 2];
        for (size_t i = 0; i < Size; ++i) {
            operand[2 * i] = ScalarType(value_ptr() + i, grad_ptr() + i);
            operand[2 * i + 1] = ScalarType(other.value_ptr() + i, other.grad_ptr() + i);
        }
        tracer.pushOperand(operand);
        return {temp, newHeadNode};
    }

    template<class PlainScalar, size_t Size>
    inline SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size>
    SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size>::operator/(const SIMD& other) const {
        auto& tracer = TracerType::getInstance();
        const auto temp = Base::operator/(other);
        const ScalarType newHeadNode = tracer.pushOperation(temp, ExprType::Div);
        ScalarType operand[Size * 2];
        for (size_t i = 0; i < Size; ++i) {
            operand[2 * i] = ScalarType(value_ptr() + i, grad_ptr() + i);
            operand[2 * i + 1] = ScalarType(other.value_ptr() + i, other.grad_ptr() + i);
        }
        tracer.pushOperand(operand);
        return {temp, newHeadNode};
    }

    template<class PlainScalar, size_t Size>
    inline SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size>
    SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size>::operator-() const {
        auto& tracer = TracerType::getInstance();
        const auto temp = Base::operator-();
        const ScalarType headNode = tracer.pushOperation(temp, ExprType::Minus);
        ScalarType operand[Size];
        for (size_t i = 0; i < Size; ++i)
            operand[i] = ScalarType(value_ptr() + i, grad_ptr() + i);
        tracer.pushOperand(operand);
        return {temp, headNode};
    }

    template<class PlainScalar, size_t Size>
    inline void SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size>::load(const ScalarType* p) {
        assert(checkContinuous(Size, p) && "[Error]: Load a uncontinuous pointer is a bug");
        headNode = *p;
        Base::load(headNode.value_ptr());
    }

    template<class PlainScalar, size_t Size>
    inline void SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size>::load_partial(int n, const ScalarType* p) {
        assert(checkContinuous(n, p) && "[Error]: Load a uncontinuous pointer is a bug");
        headNode = *p;
        Base::load_partial(n, headNode.value_ptr());
    }

    template<class PlainScalar, size_t Size>
    inline void SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size>::store(ScalarType* p) const {
        for (size_t i = 0; i < Size; ++i)
            *(p + i) = ScalarType(value_ptr() + i, grad_ptr() + i);
    }

    template<class PlainScalar, size_t Size>
    inline void SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size>::store_partial(int n, ScalarType* p) const {
        for (int i = 0; i < n; ++i)
            *(p + i) = ScalarType(value_ptr() + i, grad_ptr() + i);
    }

    template<class PlainScalar, size_t Size>
    inline typename SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size>::ScalarType
    SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size>::horizontal_add() const {
        ScalarType result = 0;
        for (size_t i = 0; i < Size; ++i)
            result += this->operator[](i);
        return result;
    }

    template<class PlainScalar, size_t Size>
    inline void SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size>::swap(SIMD& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        std::swap(headNode, obj.headNode);
    }

    template<class PlainScalar, size_t Size>
    bool SIMD<Differentiable<PlainScalar, DiffMode::Reverse, 1>, Size>::checkContinuous(int n, const ScalarType* p) {
        const PlainScalar* pValue = p->value_ptr();
        const PlainScalar* pGrad = p->grad_ptr();
        for (int i = 1; i < n; ++i) {
            if ((pValue + i) != (p + i)->value_ptr())
                return false;
            if ((pGrad + i) != (p + i)->grad_ptr())
                return false;
        }
        return true;
    }
    //////////////////////////////////////////////////////////////////
    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> mul_add(
            const SIMD<ScalarType, Size>& a,
            const SIMD<ScalarType, Size>& b,
            const SIMD<ScalarType, Size>& c) {
        if constexpr (ScalarType::isReverseDiff) {
            static_assert(Size == 2 || Size == 4 || Size == 8, "[Error]: Not implemented");
            using PlainScalar = typename ScalarType::PlainScalar;
            using TracerType = typename ScalarType::TracerType;
            auto& tracer = TracerType::getInstance();
            const auto temp = mul_add<PlainScalar, Size>(a, b, c);
            ExprType source;
            if constexpr (Size == 2)
                source = ExprType::MulAdd2;
            else if constexpr (Size == 4)
                source = ExprType::MulAdd4;
            else
                source = ExprType::MulAdd8;
            const ScalarType headNode = tracer.pushOperation(temp, source);
            tracer.pushOperand(a.getHeadNode(), b.getHeadNode(), c.getHeadNode());
            return {temp, headNode};
        }
        else
            return SIMD<ScalarType, Size>(mul_add(a.getImpl(), b.getImpl(), c.getImpl()));
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> nmul_add(
            const SIMD<ScalarType, Size> a,
            const SIMD<ScalarType, Size> b,
            const SIMD<ScalarType, Size> c) {
        static_assert(!ScalarType::isDifferentiable, "[Error]: Not implemented");
        return SIMD<ScalarType, Size>(nmul_add(a.getImpl(), b.getImpl(), c.getImpl()));
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> mul_sub(
            const SIMD<ScalarType, Size> a,
            const SIMD<ScalarType, Size> b,
            const SIMD<ScalarType, Size> c) {
        static_assert(!ScalarType::isDifferentiable, "[Error]: Not implemented");
        return SIMD<ScalarType, Size>(mul_sub(a.getImpl(), b.getImpl(), c.getImpl()));
    }
}
