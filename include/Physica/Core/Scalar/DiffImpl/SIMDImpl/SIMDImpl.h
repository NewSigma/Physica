/*
 * Copyright 2024-2026 Weibo He.
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

#include "../SIMD.h"

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order, int Size>
    SIMD<Diff<T, Mode, Order>, Size>::SIMD(int x) : values(x), grads(0) {}

    template<Scalar T, DiffMode Mode, int Order, int Size>
    SIMD<Diff<T, Mode, Order>, Size>::SIMD(double x) : values(x), grads(0) {}

    template<Scalar T, DiffMode Mode, int Order, int Size>
    SIMD<Diff<T, Mode, Order>, Size>::SIMD(Scalar auto x) {
        if constexpr (Diffable<decltype(x)>) {
            values = ValueType(x.value());
            grads = GradType(x.grad());
        }
        else {
            values = ValueType(std::move(x));
            grads = GradType(0);
        }
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    SIMD<Diff<T, Mode, Order>, Size>::SIMD(ScalarType x, int count) : values(x.value(), count), grads(x.grad()) {}

    template<Scalar T, DiffMode Mode, int Order, int Size>
    SIMD<Diff<T, Mode, Order>, Size>::SIMD(ValueType values_) : values(std::move(values_)), grads(0) {}

    template<Scalar T, DiffMode Mode, int Order, int Size>
    SIMD<Diff<T, Mode, Order>, Size>::SIMD(ValueType values_, GradType grads_) : values(std::move(values_)), grads(std::move(grads_)) {}

    template<Scalar T, DiffMode Mode, int Order, int Size>
    SIMD<Diff<T, Mode, Order>, Size>::SIMD(HalfType a, HalfType b) : values(a.value(), b.value()), grads(a.grad(), b.grad()) {}

    template<Scalar T, DiffMode Mode, int Order, int Size>
    template<int OtherOrder>
    SIMD<Diff<T, Mode, Order>, Size>::SIMD(const SIMD<Diff<T, Mode, OtherOrder>, Size>& other)
            : values(other.value()), grads(other.grad()) {}

    template<Scalar T, DiffMode Mode, int Order, int Size>
    auto SIMD<Diff<T, Mode, Order>, Size>::operator[](int index) const -> ScalarType {
        return ScalarType(values[index], grads[index]);
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    auto SIMD<Diff<T, Mode, Order>, Size>::operator+(const SIMD& other) const -> This {
        return This(values + other.values, grads + other.grads);
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    auto SIMD<Diff<T, Mode, Order>, Size>::operator-(const SIMD& other) const -> This {
        return This(values - other.values, grads - other.grads);
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    auto SIMD<Diff<T, Mode, Order>, Size>::operator*(const SIMD& x) const -> This {
        return This(values * x.value(), GradType(GradType(x) * grad() + GradType(*this) * x.grad()));
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    auto SIMD<Diff<T, Mode, Order>, Size>::operator*(const ScalarType& x) const -> This {
        return This(values * x.value(), GradType(x.template grad_mask<Order - 1>() * grad() + GradType(*this) * x.grad()));
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    auto SIMD<Diff<T, Mode, Order>, Size>::operator*(const Scalar auto& x) const -> This {
        if constexpr (Diffable<decltype(x)>)
            return *this * ScalarType(x);
        else
            return This(values * x, grads * x);
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    auto SIMD<Diff<T, Mode, Order>, Size>::operator/(const SIMD& x) const -> This {
        const auto x1 = GradType(x);
        const auto v = reciprocal(x1);
        return This(value() * v.value(), GradType((grad() * GradType(x1) - GradType(*this) * x.grad()) * square(v)));
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    auto SIMD<Diff<T, Mode, Order>, Size>::operator-() const -> This {
        return This(-values, -grads);
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    auto SIMD<Diff<T, Mode, Order>, Size>::reverse(GradType grad) const noexcept -> ValueType {
        grads += grad;
        return values;
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    void SIMD<Diff<T, Mode, Order>, Size>::load(ConstPtrTy p) noexcept {
        values.load(p.value_ptr());
        grads.load(p.grad_ptr());
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    void SIMD<Diff<T, Mode, Order>, Size>::load(ConstPtrTy p, int n) noexcept {
        values.load(p.value_ptr(), n);
        grads.load(p.grad_ptr(), n);
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    void SIMD<Diff<T, Mode, Order>, Size>::store(PtrTy p) const noexcept {
        values.store(p.value_ptr());
        grads.store(p.grad_ptr());
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    void SIMD<Diff<T, Mode, Order>, Size>::store(PtrTy p, int n) const noexcept {
        values.store(p.value_ptr(), n);
        grads.store(p.grad_ptr(), n);
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    auto SIMD<Diff<T, Mode, Order>, Size>::cutoff(int count) -> This& {
        values.cutoff(count);
        grads.cutoff(count);
        return *this;
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    auto SIMD<Diff<T, Mode, Order>, Size>::swapRealImag() const noexcept -> FullRealType {
        return FullRealType(values.swapRealImag(), grads.swapRealImag());
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    auto SIMD<Diff<T, Mode, Order>, Size>::gatherRealImag() const noexcept -> FullRealType {
        return FullRealType(values.gatherRealImag(), grads.gatherRealImag());
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    auto SIMD<Diff<T, Mode, Order>, Size>::sum() const {
        return ScalarType(values.sum(), grads.sum());
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    auto SIMD<Diff<T, Mode, Order>, Size>::max() const {
        const T x = values.max();
        return ScalarType(x, GradType::select(ValueType(x) == values, grads, GradType(0)).sum());
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    auto SIMD<Diff<T, Mode, Order>, Size>::min() const {
        const T x = values.min();
        return ScalarType(x, GradType::select(ValueType(x) == values, grads, GradType(0)).sum());
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    void SIMD<Diff<T, Mode, Order>, Size>::swap(SIMD& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        values.swap(obj.values);
        grads.swap(obj.grads);
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    auto SIMD<Diff<T, Mode, Order>, Size>::asReal() const noexcept -> FullRealType {
        if constexpr (T::isComplex)
            return FullRealType(values.asReal(), grads.asReal());
        else
            return *this;
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    auto SIMD<Diff<T, Mode, Order>, Size>::real() const noexcept -> RealType {
        return RealType(values.real(), grads.real());
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    auto SIMD<Diff<T, Mode, Order>, Size>::imag() const noexcept -> RealType {
        return RealType(values.imag(), grads.imag());
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    auto SIMD<Diff<T, Mode, Order>, Size>::zeros() noexcept -> SIMD {
        return SIMD(ValueType::zeros(), GradType::zeros());
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    auto SIMD<Diff<T, Mode, Order>, Size>::select(BoolSIMDType flags, const SIMD& x, const SIMD& y) -> This {
        return This(ValueType::select(flags, x.value(), y.value()), GradType::select(flags, x.grad(), y.grad()));
    }

    template<Scalar T, DiffMode Mode, int Order, int Size>
    auto SIMD<Diff<T, Mode, Order>, Size>::asComplex(const FullRealType& reals) -> This {
        This result{};
        result.values = ValueType::asComplex(reals.value());
        result.grads = GradType::asComplex(reals.grad());
        return result;
    }
}
