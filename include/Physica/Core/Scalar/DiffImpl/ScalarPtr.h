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
        ScalarPtr() {}
        ScalarPtr(std::pair<T*, GradPtrTy> pair_) : pair(std::move(pair_)) {}
        ScalarPtr(T* pValue, GradPtrTy pGrad) : pair(std::make_pair(pValue, pGrad)) {}
        explicit ScalarPtr(ScalarType& x);
        explicit ScalarPtr(ScalarRef<ScalarType>& x);
        ScalarPtr(const This&) = default;
        ScalarPtr(This&&) noexcept = default;
        ~ScalarPtr() = default;
        /* Operators */
        This& operator=(const This& obj) = default;
        This& operator=(This&& obj) noexcept = default;
        [[nodiscard]] inline bool operator==(const This& other) const noexcept;
        [[nodiscard]] bool operator!=(const This& other) const noexcept { return !(*this == other); }
        template<Scalar U>
        [[nodiscard]] explicit operator ScalarPtr<Diff<U, Mode, Order>>() noexcept {
            using Target = ScalarPtr<Diff<U, Mode, Order>>;
            using GradPtrTy1 = Target::GradPtrTy;
            return Target(reinterpret_cast<U*>(value_ptr()), GradPtrTy1(grad_ptr()));
        }
        [[nodiscard]] ScalarRef<ScalarType> operator*() const { return ScalarRef<ScalarType>(*this); }
        [[nodiscard]] This operator+(size_t n) { return ScalarPtr(value_ptr() + n, grad_ptr() + n); }
        inline This& operator++();
        inline This& operator--();
        inline const This operator++(int);
        inline const This operator--(int);
        [[nodiscard]] T* operator[](size_t i) const noexcept { assert(i <= Order); return arr[i]; }
        /* Operations */
        inline void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] T* value_ptr() const noexcept { return pair.first; }
        template<int GradOrder = 1>
        [[nodiscard]] auto grad_ptr() const noexcept;
    };

    template<Scalar T, DiffMode Mode, int Order>
    ScalarPtr<Diff<T, Mode, Order>>::ScalarPtr(ScalarType& x) : ScalarPtr(&x.value(), &x.grad()) {}

    template<Scalar T, DiffMode Mode, int Order>
    ScalarPtr<Diff<T, Mode, Order>>::ScalarPtr(ScalarRef<ScalarType>& x) : ScalarPtr(&x.value(), &x.grad()) {}

    template<Scalar T, DiffMode Mode, int Order>
    inline bool ScalarPtr<Diff<T, Mode, Order>>::operator==(const This& other) const noexcept {
        bool flag = pair.first == other.pair.first;
        assert(flag == (pair.second == other.pair.second) && "[Error]: Bad ScalarPtr");
        return flag;
    }

    template<Scalar T, DiffMode Mode, int Order>
    inline ScalarPtr<Diff<T, Mode, Order>>& ScalarPtr<Diff<T, Mode, Order>>::operator++() {
        for (auto& p : arr)
            p++;
        return *this;
    }

    template<Scalar T, DiffMode Mode, int Order>
    inline ScalarPtr<Diff<T, Mode, Order>>& ScalarPtr<Diff<T, Mode, Order>>::operator--() {
        for (auto& p : arr)
            p--;
        return *this;
    }

    template<Scalar T, DiffMode Mode, int Order>
    inline const ScalarPtr<Diff<T, Mode, Order>> ScalarPtr<Diff<T, Mode, Order>>::operator++(int) {
        return This(pair.first++, pair.second++);
    }

    template<Scalar T, DiffMode Mode, int Order>
    inline const ScalarPtr<Diff<T, Mode, Order>> ScalarPtr<Diff<T, Mode, Order>>::operator--(int) {
        return This(pair.first--, pair.second--);
    }

    template<Scalar T, DiffMode Mode, int Order>
    inline void ScalarPtr<Diff<T, Mode, Order>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        pair.swap(obj.pair);
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<int GradOrder>
    auto ScalarPtr<Diff<T, Mode, Order>>::grad_ptr() const noexcept {
        static_assert(GradOrder > 0, "[Error]: Invalid order");
        if constexpr (GradOrder == 1)
            return pair.second;
        else
            return pair.second.template grad_ptr<GradOrder - 1>();
    }
}
