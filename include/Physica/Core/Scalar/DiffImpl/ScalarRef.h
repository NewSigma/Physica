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

#include "../Scalar.h"

namespace Physica {
    /**
     * \class ScalarRef is a proxy type returned when dereferencing \class ScalarPtr
     */
    template<Scalar T> requires(instanceof_tx<Diff, T>)
    class ScalarRef<T> : public ScalarBase<ScalarRef<T>> {
        using This = ScalarRef<T>;
        using Base = ScalarBase<This>;

        using Tr = std::remove_const_t<T>::RealType;
        using Tv = std::remove_const_t<T>::ValueType;
        using Tg = std::remove_const_t<T>::GradType;
        constexpr static bool IsConst = std::is_const_v<T>;
        static_assert(!std::is_reference<T>::value);
    public:
        using Base::Mode;
        using Base::Order;
        using typename Base::GradType;
    private:
        ScalarPtr<T> ptr;
    public:
        ScalarRef() = default;
        __host__ __device__ explicit ScalarRef(ScalarPtr<T> ptr_) noexcept;
        ScalarRef(const This&) = default;
        ScalarRef(This&&) noexcept = default;
        ~ScalarRef() = default;
        /* Operators */
        __host__ __device__ This& operator=(const This& other);
        __host__ __device__ This& operator=(This&& other) noexcept;
        __host__ __device__ This& operator=(Scalar auto other) noexcept;
        __host__ __device__ This& operator=(int x);
        __host__ __device__ This& operator=(double x);
        __host__ __device__ void operator+=(const Scalar auto& x) noexcept;
        __host__ __device__ void operator-=(const Scalar auto& x) noexcept;
        __host__ __device__ void operator*=(const Scalar auto& x) noexcept;
        __host__ __device__ void operator/=(const Scalar auto& x) noexcept;
        [[nodiscard]] __host__ __device__ bool operator==(const This& other) const;
        using Base::operator<=>;

        [[nodiscard]] __host__ __device__ operator ScalarRef<const T>() const noexcept requires(!IsConst);
        [[nodiscard]] __host__ __device__ operator T() const requires(!ReverseDiff<T>);
        [[nodiscard]] __host__ __device__ explicit operator float() const noexcept { return float(T(*this)); }
        [[nodiscard]] __host__ __device__ explicit operator double() const noexcept { return double(T(*this)); }
        /* Operations */
        __host__ __device__ T reverse(GradType grad = 1) const noexcept;
        __host__ __device__ void zero_grad();

        __host__ __device__ void swap(This&& obj) noexcept;
        __host__ __device__ void swap(T& obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ auto value_ptr() const noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] __host__ __device__ auto grad_ptr() const noexcept;
        [[nodiscard]] __host__ __device__ auto real_ptr(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ decltype(auto) value() const noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] __host__ __device__ decltype(auto) grad() const noexcept;
    };

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ ScalarRef<T>::ScalarRef(ScalarPtr<T> ptr_) noexcept : ptr(ptr_) {}

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ auto ScalarRef<T>::operator=(const This& other) -> This& {
        value() = other.value();
        grad() = other.grad();
        return *this;
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ auto ScalarRef<T>::operator=(This&& other) noexcept -> This& {
        value() = other.value();
        grad() = other.grad();
        return *this;
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ auto ScalarRef<T>::operator=(Scalar auto other) noexcept -> This& {
        using U = std::remove_reference<decltype(other)>::type;
        Base::template static_assert_assign<U>();
        value() = other.value();
        if constexpr (Diffable<decltype(other)>)
            grad() = other.grad();
        else
            zero_grad();
        return *this;
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ auto ScalarRef<T>::operator=(int x) -> This& {
        operator=(T(x));
        return *this;
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ auto ScalarRef<T>::operator=(double x) -> This& {
        operator=(T(x));
        return *this;
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ void ScalarRef<T>::operator+=(const Scalar auto& x) noexcept {
        Base::template static_assert_assign<decltype(x)>();
        *this = *this + x;
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ void ScalarRef<T>::operator-=(const Scalar auto& x) noexcept {
        Base::template static_assert_assign<decltype(x)>();
        *this = *this - x;
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ void ScalarRef<T>::operator*=(const Scalar auto& x) noexcept {
        Base::template static_assert_assign<decltype(x)>();
        *this = *this * x;
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ void ScalarRef<T>::operator/=(const Scalar auto& x) noexcept {
        Base::template static_assert_assign<decltype(x)>();
        *this = *this / x;
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ bool ScalarRef<T>::operator==(const This& other) const {
        return value() == other.value() && grad() == other.grad();
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ ScalarRef<T>::operator ScalarRef<const T>() const noexcept requires(!IsConst) {
        return ScalarRef<const T>(ScalarPtr<T>(value_ptr(), grad_ptr()));
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ ScalarRef<T>::operator T() const requires(!ReverseDiff<T>) {
        return T(value(), grad());
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ T ScalarRef<T>::reverse(GradType grad) const noexcept {
        static_assert(Mode == DiffMode::Reverse, "[Error]: Call reverse() of a forward diff scalar is not well defined");
        decltype(auto) g = this->grad();
        const_cast<Tv&>(g.value()) += grad.value();
        if constexpr (Order != 1)
            g.reverse(grad.grad());
        return value();
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ void ScalarRef<T>::zero_grad() {
        grad() = 0;
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ void ScalarRef<T>::swap(This&& obj) noexcept {
        static_assert(!IsConst, "[Error]: Cannot modify const reference");
        std::ignore = std::move(obj); // Silent clang-tidy warning
        value().swap(obj.value());
        grad().swap(obj.grad());
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ void ScalarRef<T>::swap(T& obj) noexcept {
        static_assert(!IsConst, "[Error]: Cannot modify const reference");
        obj.swap(*this);
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ auto ScalarRef<T>::value_ptr() const noexcept {
        return ptr.value_ptr();
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    template<int GradOrder>
    __host__ __device__ auto ScalarRef<T>::grad_ptr() const noexcept {
        return ptr.template grad_ptr<GradOrder>();
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ auto ScalarRef<T>::real_ptr(this auto&& self) noexcept {
        constexpr bool ConstPtr = std::is_const_v<std::remove_reference_t<decltype(self)>>;
        using RetTy = std::conditional<ConstPtr, ScalarPtr<const Tr>, ScalarPtr<Tr>>::type;
        if constexpr (T::isComplex())
            return RetTy(self.value().real_ptr(), self.grad().real_ptr());
        else
            return RetTy(self.ptr);
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ decltype(auto) ScalarRef<T>::value() const noexcept {
        return *value_ptr();
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    template<int GradOrder>
    __host__ __device__ decltype(auto) ScalarRef<T>::grad() const noexcept {
        return *(ptr.template grad_ptr<GradOrder>());
    }

    template<Scalar T>
    std::ostream& operator<<(std::ostream& os, const ScalarRef<T>& obj) {
        return os << obj.value();
    }
}

namespace Physica {
    template<class T>
    class Traits<ScalarRef<T>> : public Traits<T> {};
}

namespace std {
    template<class T>
    struct formatter<Physica::ScalarRef<T>, char> {
        constexpr auto parse(std::format_parse_context& ctx) {
            return ctx.begin();
        }

        static auto format(const Physica::ScalarRef<T>& obj, auto& ctx) {
            return std::format_to(ctx.out(), "{}", T(obj));
        }
    };
}
