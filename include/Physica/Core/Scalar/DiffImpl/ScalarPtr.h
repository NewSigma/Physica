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
    class ScalarPtr<Diff<T, Mode, Order>> {
        using ScalarType = Diff<T, Mode, Order>;
        using This = ScalarPtr<ScalarType>;
        using RefTy = ScalarRef<ScalarType>;
        using GradType = ScalarType::GradType;
    public:
        using GradPtrTy = GradType::PtrTy;
    private:
        T* pValue;
        GradPtrTy pGrad;
    public:
        ScalarPtr() = default;
        __host__ __device__ ScalarPtr(T* pValue_, GradPtrTy pGrad_);
        __host__ __device__ explicit ScalarPtr(ScalarType& x);
        __host__ __device__ explicit ScalarPtr(RefTy& x);
        ScalarPtr(const This&) = default;
        ScalarPtr(This&&) noexcept = default;
        ~ScalarPtr() = default;
        /* Operators */
        This& operator=(const This& obj) = default;
        This& operator=(This&& obj) noexcept = default;
        [[nodiscard]] __host__ __device__ bool operator==(const This& other) const noexcept;
        [[nodiscard]] __host__ __device__ bool operator!=(const This& other) const noexcept { return !(*this == other); }
        template<Scalar U>
        [[nodiscard]] __host__ __device__ explicit operator ScalarPtr<Diff<U, Mode, Order>>() noexcept;
        [[nodiscard]] __host__ __device__ RefTy operator*() const noexcept;
        [[nodiscard]] __host__ __device__ RefTy operator->() const noexcept;
        [[nodiscard]] __host__ __device__ This operator+(size_t n) const noexcept;
        __host__ __device__ This& operator++() noexcept;
        __host__ __device__ This& operator--() noexcept;
        __host__ __device__ const This operator++(int) noexcept;
        __host__ __device__ const This operator--(int) noexcept;
        /* Operations */
        __host__ __device__ void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ T* value_ptr() const noexcept { return pValue; }
        template<int GradOrder = 1>
        [[nodiscard]] __host__ __device__ auto grad_ptr() const noexcept;
        template<size_t I>
        [[nodiscard]] __host__ __device__ constexpr T* get() const noexcept { return get(I); }
        [[nodiscard]] __host__ __device__ constexpr T* get(int order) const noexcept;
    };

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ ScalarPtr<Diff<T, Mode, Order>>::ScalarPtr(T* pValue_, GradPtrTy pGrad_) : pValue(pValue_), pGrad(pGrad_) {}

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ ScalarPtr<Diff<T, Mode, Order>>::ScalarPtr(ScalarType& x) : ScalarPtr(&x.value(), &x.grad()) {}

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ ScalarPtr<Diff<T, Mode, Order>>::ScalarPtr(RefTy& x) : ScalarPtr(&x.value(), &x.grad()) {}

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ bool ScalarPtr<Diff<T, Mode, Order>>::operator==(const This& other) const noexcept {
        bool flag = pValue == other.pValue;
        assert(flag == (pGrad == other.pGrad) && "[Error]: Bad ScalarPtr");
        return flag;
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<Scalar U>
    __host__ __device__ ScalarPtr<Diff<T, Mode, Order>>::operator ScalarPtr<Diff<U, Mode, Order>>() noexcept {
        using Target = ScalarPtr<Diff<U, Mode, Order>>;
        using GradPtrTy1 = Target::GradPtrTy;
        return Target(reinterpret_cast<U*>(value_ptr()), GradPtrTy1(grad_ptr()));
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ auto ScalarPtr<Diff<T, Mode, Order>>::operator*() const noexcept -> RefTy {
        return RefTy(*this);
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ auto ScalarPtr<Diff<T, Mode, Order>>::operator->() const noexcept -> RefTy {
        return operator*();
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ auto ScalarPtr<Diff<T, Mode, Order>>::operator+(size_t n) const noexcept -> This {
        return ScalarPtr(value_ptr() + n, grad_ptr() + n);
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ ScalarPtr<Diff<T, Mode, Order>>& ScalarPtr<Diff<T, Mode, Order>>::operator++() noexcept {
        pValue++;
        pGrad++;
        return *this;
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ ScalarPtr<Diff<T, Mode, Order>>& ScalarPtr<Diff<T, Mode, Order>>::operator--() noexcept {
        pValue--;
        pGrad--;
        return *this;
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ const ScalarPtr<Diff<T, Mode, Order>> ScalarPtr<Diff<T, Mode, Order>>::operator++(int) noexcept {
        return This(pValue++, pGrad++);
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ const ScalarPtr<Diff<T, Mode, Order>> ScalarPtr<Diff<T, Mode, Order>>::operator--(int) noexcept {
        return This(pValue--, pGrad--);
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ void ScalarPtr<Diff<T, Mode, Order>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        using std::swap;
        swap(pValue, obj.pValue);
        swap(pGrad, obj.pGrad);
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<int GradOrder>
    __host__ __device__ auto ScalarPtr<Diff<T, Mode, Order>>::grad_ptr() const noexcept {
        static_assert(GradOrder > 0, "[Error]: Invalid order");
        if constexpr (GradOrder == 1)
            return pGrad;
        else
            return pGrad.template grad_ptr<GradOrder - 1>();
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ constexpr T* ScalarPtr<Diff<T, Mode, Order>>::get(int order) const noexcept {
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

    template<Scalar T, DiffMode Mode, int Order>
    void swap(ScalarPtr<Diff<T, Mode, Order>>& p1, ScalarPtr<Diff<T, Mode, Order>>& p2) noexcept {
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
