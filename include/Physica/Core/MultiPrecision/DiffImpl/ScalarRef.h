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
    class ScalarRef<Diff<T, Mode, Order>>
            : public ScalarBase<ScalarRef<Diff<T, Mode, Order>>>
            , private ScalarPtr<Diff<T, Mode, Order>> {
        using This = ScalarRef<Diff<T, Mode, Order>>;
        using Base = ScalarBase<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::PtrTy;
        using typename Base::GradType;
    public:
        ScalarRef() = default;
        explicit ScalarRef(PtrTy ptr_) : PtrTy(ptr_) {}
        ScalarRef(const ScalarRef&) = default;
        ScalarRef(ScalarRef&&) noexcept = default;
        ~ScalarRef() = default;
        /* Operators */
        inline This& operator=(const This& other);
        inline This& operator=(const ScalarType& other);
        inline This& operator=(int x) { return operator=(ScalarType(x)); }
        inline This& operator=(double x) { return operator=(ScalarType(x)); }
        [[nodiscard]] operator ScalarType() const { return ScalarType(Base::value(), Base::grad()); }
        [[nodiscard]] __host__ __device__ explicit operator float() const noexcept { return float(ScalarType(*this)); }
        [[nodiscard]] __host__ __device__ explicit operator double() const noexcept { return double(ScalarType(*this)); }

        [[nodiscard]] ScalarType operator-() const { return -ScalarType(*this); }
        [[nodiscard]] bool operator==(const This& other) const;
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
        using PtrTy::value_ptr;
        using PtrTy::grad_ptr;
    };

    template<Scalar T, DiffMode Mode, int Order>
    inline ScalarRef<Diff<T, Mode, Order>>& ScalarRef<Diff<T, Mode, Order>>::operator=(const ScalarRef& other) {
        Base::value() = other.value();
        Base::grad() = other.grad();
        return *this;
    }

    template<Scalar T, DiffMode Mode, int Order>
    inline ScalarRef<Diff<T, Mode, Order>>& ScalarRef<Diff<T, Mode, Order>>::operator=(const ScalarType& other) {
        Base::value() = other.value();
        Base::grad() = other.grad();
        return *this;
    }

    template<Scalar T, DiffMode Mode, int Order>
    bool ScalarRef<Diff<T, Mode, Order>>::operator==(const This& other) const {
        return Base::value() == other.value() && Base::grad() == other.grad();
    }

    template<Scalar T, DiffMode Mode, int Order>
    T ScalarRef<Diff<T, Mode, Order>>::reverse(GradType grad_) const noexcept {
        static_assert(Mode == DiffMode::Reverse, "[Error]: Call reverse() of a forward diff scalar is not well defined");
        auto& g = const_cast<GradType&>(Base::grad());
        g.value() += grad_.value();
        if constexpr (Order != 1)
            g.reverse(grad_.grad());
        return Base::value();
    }

    template<Scalar T, DiffMode Mode, int Order>
    inline void ScalarRef<Diff<T, Mode, Order>>::zero_grad() {
        Base::grad() = 0;
    }

    template<Scalar T, DiffMode Mode, int Order>
    void ScalarRef<Diff<T, Mode, Order>>::swap(This&& obj) noexcept {
        Base::value().swap(obj.value());
        Base::grad().swap(obj.grad());
    }

    template<Scalar T, DiffMode Mode, int Order>
    void ScalarRef<Diff<T, Mode, Order>>::swap(ScalarType& obj) noexcept {
        obj.swap(*this);
    }

    template<Scalar T, DiffMode Mode, int Order>
    inline std::ostream& operator<<(std::ostream& os, const ScalarRef<Diff<T, Mode, Order>>& obj) {
        return os << obj.value();
    }
}
