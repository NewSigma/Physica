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
    ////////////////////////////////////////////////////////////////////////////////////
    template<class PlainScalar, size_t Size>
    SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>::SIMD(PlainScalar s) : Base(s) {
        auto& tracer = TracerType::getInstance();
        headNode = tracer.pushOperation(*this, ExprType::Set);
    }

    template<class PlainScalar, size_t Size>
    SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>::SIMD(ScalarType s) : Base(s.getValue()) {
        auto& tracer = TracerType::getInstance();
        headNode = tracer.pushOperation(*this, ExprType::Assign);
        ScalarType operand[Size];
        for (size_t i = 0; i < Size; ++i)
            operand[i] = s;
        tracer.pushOperand(operand);
    }

    template<class PlainScalar, size_t Size>
    SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>::SIMD(Base base, ScalarType headNode_)
            : Base(std::move(base)), headNode(headNode_) {}

    template<class PlainScalar, size_t Size>
    inline typename SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>::ScalarType
    SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>::operator[](int i) const {
        return ScalarType(value_ptr() + i, grad_ptr() + i);
    }

    template<class PlainScalar, size_t Size>
    inline SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>
    SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>::operator+(const SIMD& other) const {
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
    inline SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>
    SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>::operator-(const SIMD& other) const {
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
    inline SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>
    SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>::operator*(const SIMD& other) const {
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
    inline SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>
    SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>::operator*(const ScalarType& v) const {
        return operator*(SIMD(v));
    }

    template<class PlainScalar, size_t Size>
    inline SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>
    SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>::operator/(const SIMD& other) const {
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
    inline SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>
    SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>::operator-() const {
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
    inline void SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>::load(const ScalarType* p) {
        assert(checkContinuous(Size, p) && "[Error]: Load a uncontinuous pointer is a bug");
        headNode = *p;
        Base::load(headNode.value_ptr());
    }

    template<class PlainScalar, size_t Size>
    inline void SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>::load_partial(int n, const ScalarType* p) {
        assert(checkContinuous(n, p) && "[Error]: Load a uncontinuous pointer is a bug");
        headNode = *p;
        Base::load_partial(n, headNode.value_ptr());
    }

    template<class PlainScalar, size_t Size>
    inline void SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>::store(ScalarType* p) const {
        for (size_t i = 0; i < Size; ++i)
            *(p + i) = ScalarType(value_ptr() + i, grad_ptr() + i);
    }

    template<class PlainScalar, size_t Size>
    inline void SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>::store_partial(int n, ScalarType* p) const {
        for (int i = 0; i < n; ++i)
            *(p + i) = ScalarType(value_ptr() + i, grad_ptr() + i);
    }

    template<class PlainScalar, size_t Size>
    inline typename SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>::ScalarType
    SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>::sum() const {
        ScalarType result = 0;
        for (size_t i = 0; i < Size; ++i)
            result += this->operator[](i);
        return result;
    }

    template<class PlainScalar, size_t Size>
    inline void SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>::swap(SIMD& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        std::swap(headNode, obj.headNode);
    }

    template<class PlainScalar, size_t Size>
    bool SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size>::checkContinuous(int n, const ScalarType* p) {
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
}
