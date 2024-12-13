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

#include "../Diff.h"

namespace Physica::Core {
    template<Scalar T, DiffMode Mode, int Order>
    class ScalarRef<Diff<T, Mode, Order>> {
        using ScalarType = Diff<T, Mode, Order>;
        using This = ScalarRef<ScalarType>;
        using PtrTy = ScalarType::PtrTy;
        using GradType = ScalarType::GradType;
        template<int GradOrder>
        using GradRtnTy = std::conditional<Order == GradOrder, T&, ScalarRef<Diff<T, Mode, Order - GradOrder>>>::type;
    private:
        PtrTy ptr;
    public:
        explicit ScalarRef(PtrTy ptr_) : ptr(ptr_) {}
        ScalarRef(const ScalarRef&) = default;
        ScalarRef(ScalarRef&&) noexcept = default;
        ~ScalarRef() = default;
        /* Operators */
        inline ScalarRef& operator=(const This& other);
        inline ScalarRef& operator=(const ScalarType& other);
        inline ScalarRef& operator=(int x) { return operator=(ScalarType(x)); }
        inline ScalarRef& operator=(double x) { return operator=(ScalarType(x)); }
        [[nodiscard]] operator ScalarType() const { return ScalarType(getValue(), getGrad()); }
        [[nodiscard]] __host__ __device__ explicit operator float() const noexcept { return float(ScalarType(*this)); }
        [[nodiscard]] __host__ __device__ explicit operator double() const noexcept { return double(ScalarType(*this)); }
        template<Scalar U> auto operator+(const U& s) const { return ScalarType(*this) + s; }
        template<Scalar U> auto operator-(const U& s) const { return ScalarType(*this) - s; }
        template<Scalar U> auto operator*(const U& s) const { return ScalarType(*this) * s; }
        template<Scalar U> auto operator/(const U& s) const { return ScalarType(*this) / s; }
        template<Scalar U> void operator+=(const U& s) { *this = ScalarType(*this) + s; }
        template<Scalar U> void operator-=(const U& s) { *this = ScalarType(*this) - s; }
        template<Scalar U> void operator*=(const U& s) { *this = ScalarType(*this) * s; }
        template<Scalar U> void operator/=(const U& s) { *this = ScalarType(*this) / s; }
        template<Scalar U> auto operator+(const ScalarRef<U>& s) const { return operator+(U(s)); }
        template<Scalar U> auto operator-(const ScalarRef<U>& s) const { return operator-(U(s)); }
        template<Scalar U> auto operator*(const ScalarRef<U>& s) const { return operator*(U(s)); }
        template<Scalar U> auto operator/(const ScalarRef<U>& s) const { return operator/(U(s)); }
        template<Scalar U> void operator+=(const ScalarRef<U>& s) { operator+=(U(s)); }
        template<Scalar U> void operator-=(const ScalarRef<U>& s) { operator-=(U(s)); }
        template<Scalar U> void operator*=(const ScalarRef<U>& s) { operator*=(U(s)); }
        template<Scalar U> void operator/=(const ScalarRef<U>& s) { operator/=(U(s)); }
        [[nodiscard]] ScalarType operator-() const { return -ScalarType(*this); }
        __host__ __device__ inline bool operator>(double s) const noexcept { return ScalarType(*this) > s; }
        __host__ __device__ inline bool operator<(double s) const noexcept { return ScalarType(*this) < s; }
        template<Scalar U>
        __host__ __device__ bool operator>(const U& s) const noexcept { return ScalarType(*this) > s; }
        template<Scalar U>
        __host__ __device__ bool operator<(const U& s) const noexcept { return ScalarType(*this) < s; }
        template<Scalar U>
        __host__ __device__ bool operator>(const ScalarRef<U>& s) const noexcept { return operator>(U(s)); }
        template<Scalar U>
        __host__ __device__ bool operator<(const ScalarRef<U>& s) const noexcept { return operator<(U(s)); }
        /* Operations */
        T reverse(GradType grad_ = 1) const noexcept;
        inline void zero_grad();

        void swap(This&& obj) noexcept;
        void swap(ScalarType& obj) noexcept;
        /* Getters */
        [[nodiscard]] auto real() const { return ScalarType(*this).real(); }
        [[nodiscard]] auto imag() const { return ScalarType(*this).imag(); }
        [[nodiscard]] auto conjugate() const { return ScalarType(*this).conjugate(); }

        [[nodiscard]] __host__ __device__ auto* value_ptr() noexcept { return ptr.value_ptr(); }
        [[nodiscard]] __host__ __device__ const auto* value_ptr() const noexcept { return ptr.value_ptr(); }
        [[nodiscard]] __host__ __device__ auto* grad_ptr() noexcept { return ptr.grad_ptr(); }
        [[nodiscard]] __host__ __device__ const auto* grad_ptr() const noexcept { return ptr.grad_ptr(); }
        [[nodiscard]] auto& getValue() noexcept { return *ptr.value_ptr(); }
        [[nodiscard]] const auto& getValue() const noexcept { return *ptr.value_ptr(); }
        template<int GradOrder = 1>
        [[nodiscard]] GradRtnTy<GradOrder> getGrad() noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] const GradRtnTy<GradOrder> getGrad() const noexcept;
    };

    template<Scalar T, DiffMode Mode, int Order>
    inline ScalarRef<Diff<T, Mode, Order>>& ScalarRef<Diff<T, Mode, Order>>::operator=(const ScalarRef& other) {
        getValue() = other.getValue();
        getGrad() = other.getGrad();
        return *this;
    }

    template<Scalar T, DiffMode Mode, int Order>
    inline ScalarRef<Diff<T, Mode, Order>>& ScalarRef<Diff<T, Mode, Order>>::operator=(const ScalarType& other) {
        getValue() = other.getValue();
        getGrad() = other.getGrad();
        return *this;
    }

    template<Scalar T, DiffMode Mode, int Order>
    T ScalarRef<Diff<T, Mode, Order>>::reverse(GradType grad_) const noexcept {
        static_assert(Mode == DiffMode::Reverse, "[Error]: Call reverse() of a forward diff scalar is not well defined");
        auto& g = const_cast<GradType&>(getGrad());
        g.getValue() += grad_.getValue();
        if constexpr (Order != 1)
            g.reverse(grad_.getGrad());
        return getValue();
    }

    template<Scalar T, DiffMode Mode, int Order>
    inline void ScalarRef<Diff<T, Mode, Order>>::zero_grad() {
        getGrad() = 0;
    }

    template<Scalar T, DiffMode Mode, int Order>
    void ScalarRef<Diff<T, Mode, Order>>::swap(This&& obj) noexcept {
        getValue().swap(obj.getValue());
        getGrad().swap(obj.getGrad());
    }

    template<Scalar T, DiffMode Mode, int Order>
    void ScalarRef<Diff<T, Mode, Order>>::swap(ScalarType& obj) noexcept {
        obj.swap(*this);
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<int GradOrder>
    ScalarRef<Diff<T, Mode, Order>>::template GradRtnTy<GradOrder> ScalarRef<Diff<T, Mode, Order>>::getGrad() noexcept {
        if constexpr (GradOrder == 1) {
            if constexpr (Order == GradOrder)
                return *ptr.grad_ptr();
            else
                return ScalarRef<Diff<T, Mode, Order - GradOrder>>(ptr.grad_ptr());
        }
        else
            return (*ptr.grad_ptr()).template getGrad<GradOrder - 1>();
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<int GradOrder>
    const ScalarRef<Diff<T, Mode, Order>>::template GradRtnTy<GradOrder> ScalarRef<Diff<T, Mode, Order>>::getGrad() const noexcept {
        return const_cast<This&>(*this).template getGrad<GradOrder>();
    }

    template<Scalar T, DiffMode Mode, int Order>
    inline std::ostream& operator<<(std::ostream& os, const ScalarRef<Diff<T, Mode, Order>>& obj) {
        return os << obj.getValue();
    }
}
