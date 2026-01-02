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
    template<Scalar T> requires(instanceof_tx<Diff, T>)
    class ScalarPtr<T> {
        using This = ScalarPtr<T>;
        using RefTy = ScalarRef<T>;

        using Tv = std::remove_const_t<T>::ValueType;
        using Tg = std::remove_const_t<T>::GradType;
        constexpr static DiffMode Mode = ForwardDiff<T> ? DiffMode::Forward : DiffMode::Reverse;
        constexpr static int Order = Traits<T>::Order;
        constexpr static bool IsConst = std::is_const_v<T>;
    public:
        using ValuePtrTy = std::conditional<IsConst, typename Tv::ConstPtrTy, typename Tv::PtrTy>::type;
        using GradPtrTy = std::conditional<IsConst, typename Tg::ConstPtrTy, typename Tg::PtrTy>::type;
    private:
        ValuePtrTy pValue;
        GradPtrTy pGrad;
    public:
        ScalarPtr() = default;
        __host__ __device__ ScalarPtr(ValuePtrTy pValue, GradPtrTy pGrad);
        __host__ __device__ explicit ScalarPtr(T& x);
        __host__ __device__ explicit ScalarPtr(RefTy& x);
        ScalarPtr(const This&) = default;
        ScalarPtr(This&&) noexcept = default;
        ~ScalarPtr() = default;
        /* Operators */
        This& operator=(const This& obj) = default;
        This& operator=(This&& obj) noexcept = default;
        [[nodiscard]] __host__ __device__ bool operator==(const This& other) const noexcept;
        [[nodiscard]] __host__ __device__ bool operator!=(const This& other) const noexcept { return !(*this == other); }

        [[nodiscard]] __host__ __device__ operator ScalarPtr<const T>() const noexcept requires(!IsConst);
        template<Scalar U>
        [[nodiscard]] __host__ __device__ explicit operator ScalarPtr<Diff<U, Mode, Order>>() const noexcept requires(!IsConst);
        template<Scalar U>
        [[nodiscard]] __host__ __device__ explicit operator ScalarPtr<const Diff<U, Mode, Order>>() const noexcept requires(IsConst);
        [[nodiscard]] __host__ __device__ RefTy operator*() const noexcept;
        [[nodiscard]] __host__ __device__ RefTy operator->() const noexcept;
        [[nodiscard]] __host__ __device__ This operator+(size_t n) const noexcept;
        __host__ __device__ This& operator++() noexcept;
        __host__ __device__ This& operator--() noexcept;
        __host__ __device__ This operator++(int) noexcept;
        __host__ __device__ This operator--(int) noexcept;
        /* Operations */
        __host__ __device__ void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ auto value_ptr() const noexcept { return pValue; }
        template<int GradOrder = 1>
        [[nodiscard]] __host__ __device__ auto grad_ptr() const noexcept;
        template<size_t I>
        [[nodiscard]] __host__ __device__ constexpr ValuePtrTy get() const noexcept { return get(I); }
        [[nodiscard]] __host__ __device__ constexpr ValuePtrTy get(int order) const noexcept;
    };

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ ScalarPtr<T>::ScalarPtr(ValuePtrTy pValue, GradPtrTy pGrad) : pValue(pValue), pGrad(pGrad) {}

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ ScalarPtr<T>::ScalarPtr(T& x) : ScalarPtr(&x.value(), &x.grad()) {}

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ ScalarPtr<T>::ScalarPtr(RefTy& x) : ScalarPtr(&x.value(), &x.grad()) {}

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ bool ScalarPtr<T>::operator==(const This& other) const noexcept {
        bool flag = pValue == other.pValue;
        assert(flag == (pGrad == other.pGrad) && "[Error]: Bad ScalarPtr");
        return flag;
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ ScalarPtr<T>::operator ScalarPtr<const T>() const noexcept requires(!IsConst) {
        return ScalarPtr<const T>(pValue, pGrad);
    }
    /**
     * Counter part of reinterpret_cast. Use it e.g. float64 <-> cfloat64
     */
    template<Scalar T> requires(instanceof_tx<Diff, T>)
    template<Scalar U>
    __host__ __device__ ScalarPtr<T>::operator ScalarPtr<Diff<U, Mode, Order>>() const noexcept requires(!IsConst) {
        using RetTy = ScalarPtr<Diff<U, Mode, Order>>;
        using RetPtrV = RetTy::ValuePtrTy;
        using RetPtrG = RetTy::GradPtrTy;
        return RetTy(reinterpret_cast<RetPtrV>(value_ptr()), RetPtrG(grad_ptr()));
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    template<Scalar U>
    __host__ __device__ ScalarPtr<T>::operator ScalarPtr<const Diff<U, Mode, Order>>() const noexcept requires(IsConst) {
        using RetTy = ScalarPtr<const Diff<U, Mode, Order>>;
        using RetPtrV = RetTy::ValuePtrTy;
        using RetPtrG = RetTy::GradPtrTy;
        return RetTy(reinterpret_cast<RetPtrV>(value_ptr()), RetPtrG(grad_ptr()));
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ auto ScalarPtr<T>::operator*() const noexcept -> RefTy {
        return RefTy(*this);
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ auto ScalarPtr<T>::operator->() const noexcept -> RefTy {
        return operator*();
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ auto ScalarPtr<T>::operator+(size_t n) const noexcept -> This {
        return ScalarPtr(value_ptr() + n, grad_ptr() + n);
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ auto ScalarPtr<T>::operator++() noexcept -> This& {
        pValue++;
        pGrad++;
        return *this;
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ auto ScalarPtr<T>::operator--() noexcept -> This& {
        pValue--;
        pGrad--;
        return *this;
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ auto ScalarPtr<T>::operator++(int) noexcept -> This {
        return This(pValue++, pGrad++);
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ auto ScalarPtr<T>::operator--(int) noexcept -> This {
        return This(pValue--, pGrad--);
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ void ScalarPtr<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        using std::swap;
        swap(pValue, obj.pValue);
        swap(pGrad, obj.pGrad);
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    template<int GradOrder>
    __host__ __device__ auto ScalarPtr<T>::grad_ptr() const noexcept {
        static_assert(GradOrder > 0, "[Error]: Invalid order");
        if constexpr (GradOrder == 1)
            return pGrad;
        else
            return pGrad.template grad_ptr<GradOrder - 1>();
    }

    template<Scalar T> requires(instanceof_tx<Diff, T>)
    __host__ __device__ constexpr auto ScalarPtr<T>::get(int order) const noexcept -> ValuePtrTy {
        if (order == 0)
            return pValue;

        if constexpr (Order == 1) {
            assert(order == 1 && "[Error]: Invalid order");
            return pGrad;
        }
        else {
            assert(order > 0 && "[Error]: Invalid order");
            return pGrad.get(order - 1);
        }
    }

    template<Scalar T>
    void swap(ScalarPtr<T>& p1, ScalarPtr<T>& p2) noexcept {
        p1.swap(p2);
    }
}

namespace std {
    template<class T>
    struct tuple_size<Physica::ScalarPtr<T>> : public integral_constant<std::size_t, 1 + Physica::Traits<T>::Order> {};

    template<std::size_t I, class T>
    struct tuple_element<I, Physica::ScalarPtr<T>> {
        using type = Physica::Traits<T>::ValueType::PtrTy;
    };
}
