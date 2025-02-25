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

#include "../VectorExpr.h"

namespace Physica {
    template<Vector T, Scalar U>
    class VectorExpr<ExprType::Mul, T, U>
            : public BinaryVectorExpr<ExprType::Mul, T, U> {
        using Base = BinaryVectorExpr<ExprType::Mul, T, U>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isReverseDiff;
        using Base::isComplex;
    public:
        using Base::Base;
        /* Operations */
        template<Vector V, class Executor = SeqExecutor>
        inline void assign(V& v) const;

        template<Vector V, class Executor = SeqExecutor>
        inline void assign_add(V& v) const;

        [[nodiscard]] CoDiff<ScalarType> calc(size_t index) const {
            return Base::getLHS().calc(index) * Base::getRHS();
        }

        [[nodiscard]] ValueType calc_value(size_t index) const {
            return Base::getLHS().calc_value(index) * Base::getRHS().value();
        }

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const {
            return Base::getLHS().template packet<Pack>(index) * Pack(Base::getRHS());
        }

        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const {
            return Base::getLHS().template packetPartial<Pack>(index, count) * Pack(Base::getRHS());
        }

        [[nodiscard]] ScalarType sum() const { return Base::getLHS().sum() * Base::getRHS(); }

        template<Vector V>
        void reverse(const V& grad_) const noexcept requires(isReverseDiff) {
            const auto& grad = grad_.values();
            const auto& lhs = Base::getLHS();
            const auto& rhs = Base::getRHS();
            if constexpr (ReverseDiff<T>)
                lhs.reverse(rhs.value() * grad);
            if constexpr (ReverseDiff<U>)
                rhs.reverse(lhs.values() * grad);
        }
    private:
        void assign_add_mkl(void* y) const noexcept;
    };

    template<Vector T, Scalar U>
    template<Vector V, class Executor>
    inline void VectorExpr<ExprType::Mul, T, U>::assign(V& v) const {
        constexpr bool FastAssign = Traits<std::remove_cvref_t<T>>::FastAssign;
        if constexpr (FastAssign) {
            Base::getLHS().template assign<V, Executor>(v);
            v *= Base::getRHS();
        }
        else
            Base::assign(v);
    }

    template<Vector T, Scalar U>
    template<Vector V, class Executor>
    inline void VectorExpr<ExprType::Mul, T, U>::assign_add(V& v) const {
        if constexpr (Internal::EnableMKL<T, V>::value)
            assign_add_mkl(v.data());
        else
            Base::assign_add(v);
    }

    template<Vector T1, Vector T2>
    class VectorExpr<ExprType::Mul, T1, T2>
            : public BinaryVectorExpr<ExprType::Mul, T1, T2> {
        using Base = BinaryVectorExpr<ExprType::Mul, T1, T2>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] CoDiff<ScalarType> calc(size_t index) const {
            return Base::getLHS().calc(index) * Base::getRHS().calc(index);
        }

        [[nodiscard]] ValueType calc_value(size_t index) const {
            return Base::getLHS().calc_value(index) * Base::getRHS().calc_value(index);
        }

        template<Packet Pack>
        [[nodiscard]] Pack packet(size_t index) const {
            return Base::getLHS().template packet<Pack>(index) * Base::getRHS().template packet<Pack>(index);
        }

        template<Packet Pack>
        [[nodiscard]] Pack packetPartial(size_t index, size_t count) const {
            return Base::getLHS().template packetPartial<Pack>(index, count) * Base::getRHS().template packetPartial<Pack>(index, count);
        }

        template<class U>
        void reverse(const U& grad_) const noexcept requires(isReverseDiff);
    };

    template<Vector T1, Vector T2>
    template<class U>
    void VectorExpr<ExprType::Mul, T1, T2>::reverse(const U& grad_) const noexcept requires(isReverseDiff) {
        if constexpr (Scalar<U>) {
            if constexpr (ReverseDiff<T1>)
                Base::getLHS().reverse(Base::getRHS().values() * grad_);
            if constexpr (ReverseDiff<T2>)
                Base::getRHS().reverse(Base::getLHS().values() * grad_);
        }
        else {
            static_assert(Vector<U>);
            const auto& grad = grad_.values();
            if constexpr (ReverseDiff<T1>)
                Base::getLHS().reverse(hadamard(Base::getRHS().values(), grad));
            if constexpr (ReverseDiff<T2>)
                Base::getRHS().reverse(hadamard(Base::getLHS().values(), grad));
        }
    }

    template<Vector T, Scalar U>
    [[nodiscard]] inline auto operator*(T&& v, U&& x) noexcept requires(!CUDA<T>) {
        return VectorExpr<ExprType::Mul, T&&, U&&>(std::forward<T>(v), std::forward<U>(x));
    }

    template<Scalar U, Vector T>
    [[nodiscard]] inline auto operator*(U&& x, T&& v) noexcept requires(!CUDA<T>) {
        return std::forward<T>(v) * std::forward<U>(x);
    }

    template<Vector T1, Vector T2>
    [[nodiscard]] inline auto hadamard(T1&& v1, T2&& v2) noexcept requires(!CUDA<T1> && !CUDA<T2>) {
        return VectorExpr<ExprType::Mul, T1&&, T2&&>(std::forward<T1>(v1), std::forward<T2>(v2));
    }
}

#include "MKL/Mul.h"
