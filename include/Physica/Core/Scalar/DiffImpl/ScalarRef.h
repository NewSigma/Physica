/*
 * Copyright 2024-2025 Weibo He.
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

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order>
    class ScalarRef<Diff<T, Mode, Order>>
            : public ScalarBase<ScalarRef<Diff<T, Mode, Order>>>
            , private ScalarPtr<Diff<T, Mode, Order>> {
        using This = ScalarRef<Diff<T, Mode, Order>>;
        using Base = ScalarBase<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using typename Base::PtrTy;
        using typename Base::GradType;
    private:
        template<int GradOrder>
        using GradRtnTy = std::conditional<Order == GradOrder, ValueType&, ScalarRef<Diff<ValueType, Mode, Order - GradOrder>>>::type;
    public:
        ScalarRef() = default;
        __host__ __device__ explicit ScalarRef(PtrTy ptr_) : PtrTy(ptr_) {}
        ScalarRef(const This&) = default;
        ScalarRef(This&&) noexcept = default;
        ~ScalarRef() = default;
        /* Operators */
        __host__ __device__ inline This& operator=(const This& other);
        __host__ __device__ inline This& operator=(This&& other) noexcept;
        template<Scalar U>
        __host__ __device__ inline This& operator=(const U& other);
        __host__ __device__ inline This& operator=(int x) { return operator=(T(x)); }
        __host__ __device__ inline This& operator=(double x) { return operator=(T(x)); }
        [[nodiscard]] __host__ __device__ operator ScalarType() const requires(!ReverseDiff<T>);
        [[nodiscard]] __host__ __device__ explicit operator float() const noexcept { return float(ScalarType(*this)); }
        [[nodiscard]] __host__ __device__ explicit operator double() const noexcept { return double(ScalarType(*this)); }

        [[nodiscard]] __host__ __device__ ScalarType operator-() const { return -ScalarType(*this); }
        [[nodiscard]] __host__ __device__ bool operator==(const This& other) const;
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
        __host__ __device__ T reverse(GradType grad_ = 1) const noexcept;
        __host__ __device__ inline void zero_grad();

        __host__ __device__ void swap(This&& obj) noexcept;
        __host__ __device__ void swap(ScalarType& obj) noexcept;
        /* Getters */
        using PtrTy::value_ptr;
        using PtrTy::grad_ptr;
        using Base::value;
        template<int GradOrder = 1>
        [[nodiscard]] __host__ __device__ GradRtnTy<GradOrder> grad() noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] __host__ __device__ const GradRtnTy<GradOrder> grad() const noexcept;
    };

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ inline auto ScalarRef<Diff<T, Mode, Order>>::operator=(const This& other) -> This& {
        value() = other.value();
        grad() = other.grad();
        return *this;
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ inline auto ScalarRef<Diff<T, Mode, Order>>::operator=(This&& other) noexcept -> This& {
        value() = other.value();
        grad() = other.grad();
        return *this;
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<Scalar U>
    __host__ __device__ inline auto ScalarRef<Diff<T, Mode, Order>>::operator=(const U& other) -> This& {
        value() = other.value();
        if constexpr (Diffable<U>)
            grad() = other.grad();
        else
            zero_grad();
        return *this;
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ ScalarRef<Diff<T, Mode, Order>>::operator ScalarType() const requires(!ReverseDiff<T>) {
        return ScalarType(value(), grad());
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ bool ScalarRef<Diff<T, Mode, Order>>::operator==(const This& other) const {
        return value() == other.value() && grad() == other.grad();
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ T ScalarRef<Diff<T, Mode, Order>>::reverse(GradType grad_) const noexcept {
        static_assert(Mode == DiffMode::Reverse, "[Error]: Call reverse() of a forward diff scalar is not well defined");
        auto& g = const_cast<GradType&>(grad());
        g.value() += grad_.value();
        if constexpr (Order != 1)
            g.reverse(grad_.grad());
        return value();
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ inline void ScalarRef<Diff<T, Mode, Order>>::zero_grad() {
        grad() = 0;
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ void ScalarRef<Diff<T, Mode, Order>>::swap(This&& obj) noexcept {
        value().swap(obj.value());
        grad().swap(obj.grad());
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ void ScalarRef<Diff<T, Mode, Order>>::swap(ScalarType& obj) noexcept {
        obj.swap(*this);
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<int GradOrder>
    __host__ __device__ auto ScalarRef<Diff<T, Mode, Order>>::grad() noexcept -> GradRtnTy<GradOrder> {
        return *grad_ptr();
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<int GradOrder>
    __host__ __device__ auto ScalarRef<Diff<T, Mode, Order>>::grad() const noexcept -> const GradRtnTy<GradOrder> {
        return *grad_ptr();
    }

    template<Scalar T, DiffMode Mode, int Order>
    inline std::ostream& operator<<(std::ostream& os, const ScalarRef<Diff<T, Mode, Order>>& obj) {
        return os << obj.value();
    }
}

namespace std {
    template<class T>
    struct formatter<Physica::ScalarRef<T>, char> {
        constexpr auto parse(std::format_parse_context& ctx) {
            return formatter<T, char>::parse(ctx);
        }

        auto format(const Physica::ScalarRef<T>& obj, std::format_context& ctx) const {
            return formatter<T, char>::format(T(obj), ctx);
        }
    };
}
