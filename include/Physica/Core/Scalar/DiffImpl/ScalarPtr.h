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
        using GradType = ScalarType::GradType;
    public:
        using GradPtrTy = GradType::PtrTy;
    private:
        union {
            std::pair<T*, GradPtrTy> pair;
            T* arr[Order + 1];
        };
    public:
        __host__ __device__ ScalarPtr() {}
        __host__ __device__ ScalarPtr(std::pair<T*, GradPtrTy> pair_) : pair(std::move(pair_)) {}
        __host__ __device__ ScalarPtr(T* pValue, GradPtrTy pGrad) : pair(std::make_pair(pValue, pGrad)) {}
        __host__ __device__ explicit ScalarPtr(ScalarType& x);
        __host__ __device__ explicit ScalarPtr(ScalarRef<ScalarType>& x);
        ScalarPtr(const This&) = default;
        ScalarPtr(This&&) noexcept = default;
        ~ScalarPtr() = default;
        /* Operators */
        This& operator=(const This& obj) = default;
        This& operator=(This&& obj) noexcept = default;
        [[nodiscard]] __host__ __device__ inline bool operator==(const This& other) const noexcept;
        [[nodiscard]] __host__ __device__ bool operator!=(const This& other) const noexcept { return !(*this == other); }
        template<Scalar U>
        [[nodiscard]] __host__ __device__ explicit operator ScalarPtr<Diff<U, Mode, Order>>() noexcept;
        [[nodiscard]] __host__ __device__ ScalarRef<ScalarType> operator*() const;
        [[nodiscard]] __host__ __device__ This operator+(size_t n);
        __host__ __device__ inline This& operator++();
        __host__ __device__ inline This& operator--();
        __host__ __device__ inline const This operator++(int);
        __host__ __device__ inline const This operator--(int);
        [[nodiscard]] __host__ __device__ T* operator[](size_t i) const noexcept;
        /* Operations */
        __host__ __device__ inline void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ T* value_ptr() const noexcept { return pair.first; }
        template<int GradOrder = 1>
        [[nodiscard]] __host__ __device__ auto grad_ptr() const noexcept;
    };

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ ScalarPtr<Diff<T, Mode, Order>>::ScalarPtr(ScalarType& x) : ScalarPtr(&x.value(), &x.grad()) {}

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ ScalarPtr<Diff<T, Mode, Order>>::ScalarPtr(ScalarRef<ScalarType>& x) : ScalarPtr(&x.value(), &x.grad()) {}

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ inline bool ScalarPtr<Diff<T, Mode, Order>>::operator==(const This& other) const noexcept {
        bool flag = pair.first == other.pair.first;
        assert(flag == (pair.second == other.pair.second) && "[Error]: Bad ScalarPtr");
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
    __host__ __device__ auto ScalarPtr<Diff<T, Mode, Order>>::operator*() const -> ScalarRef<ScalarType> {
        return ScalarRef<ScalarType>(*this);
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ auto ScalarPtr<Diff<T, Mode, Order>>::operator+(size_t n) -> This {
        return ScalarPtr(value_ptr() + n, grad_ptr() + n);
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ inline ScalarPtr<Diff<T, Mode, Order>>& ScalarPtr<Diff<T, Mode, Order>>::operator++() {
        for (auto& p : arr)
            p++;
        return *this;
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ inline ScalarPtr<Diff<T, Mode, Order>>& ScalarPtr<Diff<T, Mode, Order>>::operator--() {
        for (auto& p : arr)
            p--;
        return *this;
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ inline const ScalarPtr<Diff<T, Mode, Order>> ScalarPtr<Diff<T, Mode, Order>>::operator++(int) {
        return This(pair.first++, pair.second++);
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ inline const ScalarPtr<Diff<T, Mode, Order>> ScalarPtr<Diff<T, Mode, Order>>::operator--(int) {
        return This(pair.first--, pair.second--);
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ T* ScalarPtr<Diff<T, Mode, Order>>::operator[](size_t i) const noexcept {
        assert(i <= Order);
        return arr[i];
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ inline void ScalarPtr<Diff<T, Mode, Order>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        pair.swap(obj.pair);
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<int GradOrder>
    __host__ __device__ auto ScalarPtr<Diff<T, Mode, Order>>::grad_ptr() const noexcept {
        static_assert(GradOrder > 0, "[Error]: Invalid order");
        if constexpr (GradOrder == 1)
            return pair.second;
        else
            return pair.second.template grad_ptr<GradOrder - 1>();
    }
}
