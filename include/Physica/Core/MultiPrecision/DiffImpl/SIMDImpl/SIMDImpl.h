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
    template<class T, int Order, size_t Size>
    SIMD<Diff<T, DiffMode::Forward, Order>, Size>::SIMD(int x) : values(x), grads(0) {}

    template<class T, int Order, size_t Size>
    SIMD<Diff<T, DiffMode::Forward, Order>, Size>::SIMD(double x) : values(x), grads(0) {}

    template<class T, int Order, size_t Size>
    SIMD<Diff<T, DiffMode::Forward, Order>, Size>::SIMD(ScalarType x) : values(x.getValue()), grads(x.getGrad()) {}

    template<class T, int Order, size_t Size>
    SIMD<Diff<T, DiffMode::Forward, Order>, Size>::SIMD(ScalarType x, int count) : values(x.getValue(), count), grads(x.getGrad()) {}

    template<class T, int Order, size_t Size>
    SIMD<Diff<T, DiffMode::Forward, Order>, Size>::SIMD(ValuePacket values_) : values(std::move(values_)), grads(0) {}

    template<class T, int Order, size_t Size>
    SIMD<Diff<T, DiffMode::Forward, Order>, Size>::SIMD(ValuePacket values_, GradPacket grads_) : values(std::move(values_)), grads(std::move(grads_)) {}

    template<class T, int Order, size_t Size>
    SIMD<Diff<T, DiffMode::Forward, Order>, Size>::SIMD(HalfType a, HalfType b) : values(a.getValue(), b.getValue()), grads(a.getGrad(), b.getGrad()) {}

    template<class T, int Order, size_t Size>
    template<int OtherOrder>
    SIMD<Diff<T, DiffMode::Forward, Order>, Size>::SIMD(const SIMD<Diff<T, DiffMode::Forward, OtherOrder>, Size>& other)
            : values(other.getValue()), grads(other.getGrad()) {}

    template<class T, int Order, size_t Size>
    inline typename SIMD<Diff<T, DiffMode::Forward, Order>, Size>::ScalarType
    SIMD<Diff<T, DiffMode::Forward, Order>, Size>::operator[](int index) const {
        return ScalarType(values[index], grads[index]);
    }

    template<class T, int Order, size_t Size>
    inline SIMD<Diff<T, DiffMode::Forward, Order>, Size> SIMD<Diff<T, DiffMode::Forward, Order>, Size>::operator+(const SIMD& other) const {
        return This(values + other.values, grads + other.grads);
    }

    template<class T, int Order, size_t Size>
    inline SIMD<Diff<T, DiffMode::Forward, Order>, Size> SIMD<Diff<T, DiffMode::Forward, Order>, Size>::operator-(const SIMD& other) const {
        return This(values - other.values, grads - other.grads);
    }

    template<class T, int Order, size_t Size>
    inline SIMD<Diff<T, DiffMode::Forward, Order>, Size> SIMD<Diff<T, DiffMode::Forward, Order>, Size>::operator*(const SIMD& x) const {
        return This(values * x.getValue(), GradPacket(GradPacket(x) * getGrad() + GradPacket(*this) * x.getGrad()));
    }

    template<class T, int Order, size_t Size>
    inline SIMD<Diff<T, DiffMode::Forward, Order>, Size> SIMD<Diff<T, DiffMode::Forward, Order>, Size>::operator*(const ScalarType& x) const {
        return This(values * x.getValue(), GradPacket(GradType(x) * getGrad() + GradPacket(*this) * x.getGrad()));
    }

    template<class T, int Order, size_t Size>
    template<class OtherScalar>
    inline SIMD<Diff<T, DiffMode::Forward, Order>, Size> SIMD<Diff<T, DiffMode::Forward, Order>, Size>::operator*(const ScalarBase<OtherScalar>& x_) const {
        const auto& x = x_.getDerived();
        if constexpr (OtherScalar::isDifferentiable)
            return *this * ScalarType(x);
        else
            return This(values * x, grads * x);
    }

    template<class T, int Order, size_t Size>
    inline SIMD<Diff<T, DiffMode::Forward, Order>, Size> SIMD<Diff<T, DiffMode::Forward, Order>, Size>::operator/(const SIMD& x) const {
        const auto x1 = GradPacket(x);
        const auto v = reciprocal(x1);
        return This(getValue() * v.getValue(), GradPacket((getGrad() * GradPacket(x1) - GradPacket(*this) * x.getGrad()) * square(v)));
    }

    template<class T, int Order, size_t Size>
    inline SIMD<Diff<T, DiffMode::Forward, Order>, Size> SIMD<Diff<T, DiffMode::Forward, Order>, Size>::operator-() const {
        return This(-values, -grads);
    }

    template<class T, int Order, size_t Size>
    inline void SIMD<Diff<T, DiffMode::Forward, Order>, Size>::load(ConstPtrTy p) {
        values.load(p.value_ptr());
        grads.load(p.grad_ptr());
    }

    template<class T, int Order, size_t Size>
    inline void SIMD<Diff<T, DiffMode::Forward, Order>, Size>::load_partial(ConstPtrTy p, int n) {
        values.load_partial(p.value_ptr(), n);
        grads.load_partial(p.grad_ptr(), n);
    }

    template<class T, int Order, size_t Size>
    inline void SIMD<Diff<T, DiffMode::Forward, Order>, Size>::store(PtrTy p) const {
        values.store(p.value_ptr());
        grads.store(p.grad_ptr());
    }

    template<class T, int Order, size_t Size>
    inline void SIMD<Diff<T, DiffMode::Forward, Order>, Size>::store_partial(PtrTy p, int n) const {
        values.store_partial(p.value_ptr(), n);
        grads.store_partial(p.grad_ptr(), n);
    }

    template<class T, int Order, size_t Size>
    inline SIMD<Diff<T, DiffMode::Forward, Order>, Size>& SIMD<Diff<T, DiffMode::Forward, Order>, Size>::cutoff(int count) {
        values.cutoff(count);
        grads.cutoff(count);
        return *this;
    }

    template<class T, int Order, size_t Size>
    inline typename SIMD<Diff<T, DiffMode::Forward, Order>, Size>::AsRealRtnTy
    SIMD<Diff<T, DiffMode::Forward, Order>, Size>::swapRealImag() const noexcept {
        return AsRealRtnTy(values.swapRealImag(), grads.swapRealImag());
    }

    template<class T, int Order, size_t Size>
    inline typename SIMD<Diff<T, DiffMode::Forward, Order>, Size>::AsRealRtnTy
    SIMD<Diff<T, DiffMode::Forward, Order>, Size>::permRealImag() const noexcept {
        return AsRealRtnTy(values.permRealImag(), grads.permRealImag());
    }

    template<class T, int Order, size_t Size>
    inline Diff<T, DiffMode::Forward, Order> SIMD<Diff<T, DiffMode::Forward, Order>, Size>::sum() const {
        return ScalarType(values.sum(), grads.sum());
    }

    template<class T, int Order, size_t Size>
    inline Diff<T, DiffMode::Forward, Order> SIMD<Diff<T, DiffMode::Forward, Order>, Size>::max() const {
        const T x = values.max();
        return ScalarType(x, GradPacket::select(ValuePacket(x) == values, grads, GradPacket(0)).sum());
    }

    template<class T, int Order, size_t Size>
    inline Diff<T, DiffMode::Forward, Order> SIMD<Diff<T, DiffMode::Forward, Order>, Size>::min() const {
        const T x = values.min();
        return ScalarType(x, GradPacket::select(ValuePacket(x) == values, grads, GradPacket(0)).sum());
    }

    template<class T, int Order, size_t Size>
    void SIMD<Diff<T, DiffMode::Forward, Order>, Size>::swap(SIMD& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        values.swap(obj.values);
        grads.swap(obj.grads);
    }

    template<class T, int Order, size_t Size>
    typename SIMD<Diff<T, DiffMode::Forward, Order>, Size>::AsRealRtnTy SIMD<Diff<T, DiffMode::Forward, Order>, Size>::asReal() const noexcept {
        if constexpr (T::isComplex)
            return AsRealRtnTy(values.asReal(), grads.asReal());
        else
            return *this;
    }

    template<class T, int Order, size_t Size>
    inline typename SIMD<Diff<T, DiffMode::Forward, Order>, Size>::RealType
    SIMD<Diff<T, DiffMode::Forward, Order>, Size>::real() const noexcept {
        return RealType(values.real(), grads.real());
    }

    template<class T, int Order, size_t Size>
    inline typename SIMD<Diff<T, DiffMode::Forward, Order>, Size>::RealType
    SIMD<Diff<T, DiffMode::Forward, Order>, Size>::imag() const noexcept {
        return RealType(values.imag(), grads.imag());
    }

    template<class T, int Order, size_t Size>
    SIMD<Diff<T, DiffMode::Forward, Order>, Size> SIMD<Diff<T, DiffMode::Forward, Order>, Size>::select(BoolSIMDType flags, const SIMD& x, const SIMD& y) {
        return This(ValuePacket::select(flags, x.getValue(), y.getValue()), GradPacket::select(flags, x.getGrad(), y.getGrad()));
    }

    template<class T, int Order, size_t Size>
    SIMD<Diff<T, DiffMode::Forward, Order>, Size> SIMD<Diff<T, DiffMode::Forward, Order>, Size>::asComplex(const AsRealRtnTy& reals) {
        This result{};
        result.values = ValuePacket::asComplex(reals.getValue());
        result.grads = GradPacket::asComplex(reals.getGrad());
        return result;
    }
    ////////////////////////////////////////////////////////////////////////////////////
    template<class T, size_t Size>
    SIMD<Diff<T, DiffMode::Reverse, 1>, Size>::SIMD(T s) : Base(s) {
        auto& tracer = TracerType::getInstance();
        headNode = tracer.pushOperation(*this, ExprType::Set);
    }

    template<class T, size_t Size>
    SIMD<Diff<T, DiffMode::Reverse, 1>, Size>::SIMD(ScalarType s) : Base(s.getValue()) {
        auto& tracer = TracerType::getInstance();
        headNode = tracer.pushOperation(*this, ExprType::Assign);
        ScalarType operand[Size];
        for (size_t i = 0; i < Size; ++i)
            operand[i] = s;
        tracer.pushOperand(operand);
    }

    template<class T, size_t Size>
    SIMD<Diff<T, DiffMode::Reverse, 1>, Size>::SIMD(Base base, ScalarType headNode_)
            : Base(std::move(base)), headNode(headNode_) {}

    template<class T, size_t Size>
    inline typename SIMD<Diff<T, DiffMode::Reverse, 1>, Size>::ScalarType
    SIMD<Diff<T, DiffMode::Reverse, 1>, Size>::operator[](int i) const {
        return ScalarType(value_ptr() + i, grad_ptr() + i);
    }

    template<class T, size_t Size>
    inline SIMD<Diff<T, DiffMode::Reverse, 1>, Size>
    SIMD<Diff<T, DiffMode::Reverse, 1>, Size>::operator+(const SIMD& other) const {
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

    template<class T, size_t Size>
    inline SIMD<Diff<T, DiffMode::Reverse, 1>, Size>
    SIMD<Diff<T, DiffMode::Reverse, 1>, Size>::operator-(const SIMD& other) const {
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

    template<class T, size_t Size>
    inline SIMD<Diff<T, DiffMode::Reverse, 1>, Size>
    SIMD<Diff<T, DiffMode::Reverse, 1>, Size>::operator*(const SIMD& other) const {
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

    template<class T, size_t Size>
    inline SIMD<Diff<T, DiffMode::Reverse, 1>, Size>
    SIMD<Diff<T, DiffMode::Reverse, 1>, Size>::operator*(const ScalarType& v) const {
        return operator*(SIMD(v));
    }

    template<class T, size_t Size>
    inline SIMD<Diff<T, DiffMode::Reverse, 1>, Size>
    SIMD<Diff<T, DiffMode::Reverse, 1>, Size>::operator/(const SIMD& other) const {
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

    template<class T, size_t Size>
    inline SIMD<Diff<T, DiffMode::Reverse, 1>, Size>
    SIMD<Diff<T, DiffMode::Reverse, 1>, Size>::operator-() const {
        auto& tracer = TracerType::getInstance();
        const auto temp = Base::operator-();
        const ScalarType headNode = tracer.pushOperation(temp, ExprType::Minus);
        ScalarType operand[Size];
        for (size_t i = 0; i < Size; ++i)
            operand[i] = ScalarType(value_ptr() + i, grad_ptr() + i);
        tracer.pushOperand(operand);
        return {temp, headNode};
    }

    template<class T, size_t Size>
    inline void SIMD<Diff<T, DiffMode::Reverse, 1>, Size>::load(const ScalarType* p) {
        assert(checkContinuous(p, Size) && "[Error]: Load a uncontinuous pointer is a bug");
        headNode = *p;
        Base::load(headNode.value_ptr());
    }

    template<class T, size_t Size>
    inline void SIMD<Diff<T, DiffMode::Reverse, 1>, Size>::load_partial(const ScalarType* p, int n) {
        assert(checkContinuous(p, n) && "[Error]: Load a uncontinuous pointer is a bug");
        headNode = *p;
        Base::load_partial(headNode.value_ptr(), n);
    }

    template<class T, size_t Size>
    inline void SIMD<Diff<T, DiffMode::Reverse, 1>, Size>::store(ScalarType* p) const {
        for (size_t i = 0; i < Size; ++i)
            *(p + i) = ScalarType(value_ptr() + i, grad_ptr() + i);
    }

    template<class T, size_t Size>
    inline void SIMD<Diff<T, DiffMode::Reverse, 1>, Size>::store_partial(ScalarType* p, int n) const {
        for (int i = 0; i < n; ++i)
            *(p + i) = ScalarType(value_ptr() + i, grad_ptr() + i);
    }

    template<class T, size_t Size>
    inline typename SIMD<Diff<T, DiffMode::Reverse, 1>, Size>::ScalarType
    SIMD<Diff<T, DiffMode::Reverse, 1>, Size>::sum() const {
        ScalarType result = 0;
        for (size_t i = 0; i < Size; ++i)
            result += this->operator[](i);
        return result;
    }

    template<class T, size_t Size>
    inline void SIMD<Diff<T, DiffMode::Reverse, 1>, Size>::swap(SIMD& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        std::swap(headNode, obj.headNode);
    }

    template<class T, size_t Size>
    bool SIMD<Diff<T, DiffMode::Reverse, 1>, Size>::checkContinuous(const ScalarType* p, int n) {
        const T* pValue = p->value_ptr();
        const T* pGrad = p->grad_ptr();
        for (int i = 1; i < n; ++i) {
            if ((pValue + i) != (p + i)->value_ptr())
                return false;
            if ((pGrad + i) != (p + i)->grad_ptr())
                return false;
        }
        return true;
    }
    ////////////////////////////////////////////////////////////////////////////////////
    template<class T, DiffMode Mode, int Order, size_t Size>
    [[nodiscard]] inline SIMD<Diff<T, Mode, Order>, Size> mul_add(
            const SIMD<Diff<T, Mode, Order>, Size>& a,
            const SIMD<Diff<T, Mode, Order>, Size>& b,
            const SIMD<Diff<T, Mode, Order>, Size>& c) {
        using ResultType = SIMD<Diff<T, Mode, Order>, Size>;
        if constexpr (Mode == DiffMode::Forward) {
            using GradPacket = typename ResultType::GradPacket;
            auto value = mul_add(a.getValue(), b.getValue(), c.getValue());
            auto grad1 = mul_add(GradPacket(a), b.getGrad(), c.getGrad());
            auto grad2 = mul_add(GradPacket(b), a.getGrad(), grad1);
            return ResultType(std::move(value), std::move(grad2));
        }
        else {
            static_assert(Size == 2 || Size == 4 || Size == 8, "[Error]: Not implemented");
            using TracerType = typename ResultType::TracerType;
            auto& tracer = TracerType::getInstance();
            const auto temp = mul_add<T, Size>(a, b, c);
            ExprType source;
            if constexpr (Size == 2)
                source = ExprType::MulAdd2;
            else if constexpr (Size == 4)
                source = ExprType::MulAdd4;
            else
                source = ExprType::MulAdd8;
            const auto headNode = tracer.pushOperation(temp, source);
            tracer.pushOperand(a.getHeadNode(), b.getHeadNode(), c.getHeadNode());
            return {temp, headNode};
        }
    }
}
