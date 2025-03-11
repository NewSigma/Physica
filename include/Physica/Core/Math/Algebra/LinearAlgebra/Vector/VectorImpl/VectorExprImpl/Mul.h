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
    template<Vector V, Scalar U>
    class VectorExpr<ExprType::Mul, V, U>
            : public BinaryVectorExpr<ExprType::Mul, V, U> {
        using Base = BinaryVectorExpr<ExprType::Mul, V, U>;
    public:
        using Base::isReverseDiff;
        using Base::isComplex;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        template<Vector V1, class Executor = SeqExecutor>
        inline void assign(V1& v) const;

        template<Vector V1, class Executor = SeqExecutor>
        inline void assign_add(V1& v) const;

        [[nodiscard]] CoDiff<T> calc(size_t index) const {
            return Base::getLHS().calc(index) * Base::getRHS();
        }

        [[nodiscard]] Tv calc_value(size_t index) const {
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

        [[nodiscard]] T sum() const { return Base::getLHS().sum() * Base::getRHS(); }

        template<Vector V1>
        void reverse(const V1& grad_) const noexcept requires(isReverseDiff) {
            const auto& grad = grad_.values();
            const auto& lhs = Base::getLHS();
            const auto& rhs = Base::getRHS();
            if constexpr (ReverseDiff<V>)
                lhs.reverse(rhs.value() * grad);
            if constexpr (ReverseDiff<U>)
                rhs.reverse(lhs.values() * grad);
        }
    private:
        void assign_add_mkl(void* y) const noexcept;
    };

    template<Vector V, Scalar U>
    template<Vector V1, class Executor>
    inline void VectorExpr<ExprType::Mul, V, U>::assign(V1& v) const {
        constexpr bool FastAssign = Traits<std::remove_cvref_t<V>>::FastAssign;
        if constexpr (FastAssign) {
            Base::getLHS().template assign<V1, Executor>(v);
            v *= Base::getRHS();
        }
        else
            Base::assign(v);
    }

    template<Vector V, Scalar U>
    template<Vector V1, class Executor>
    inline void VectorExpr<ExprType::Mul, V, U>::assign_add(V1& v) const {
        if constexpr (Internal::EnableMKL<V, V1>::value)
            assign_add_mkl(v.data());
        else
            Base::assign_add(v);
    }

    template<Vector V1, Vector V2>
    class VectorExpr<ExprType::Mul, V1, V2>
            : public BinaryVectorExpr<ExprType::Mul, V1, V2> {
        using Base = BinaryVectorExpr<ExprType::Mul, V1, V2>;
    public:
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    public:
        using Base::Base;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t index) const {
            return Base::getLHS().calc(index) * Base::getRHS().calc(index);
        }

        [[nodiscard]] Tv calc_value(size_t index) const {
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

    template<Vector V1, Vector V2>
    template<class U>
    void VectorExpr<ExprType::Mul, V1, V2>::reverse(const U& grad_) const noexcept requires(isReverseDiff) {
        if constexpr (Scalar<U>) {
            if constexpr (ReverseDiff<V1>)
                Base::getLHS().reverse(Base::getRHS().values() * grad_);
            if constexpr (ReverseDiff<V2>)
                Base::getRHS().reverse(Base::getLHS().values() * grad_);
        }
        else {
            static_assert(Vector<U>);
            const auto& grad = grad_.values();
            if constexpr (ReverseDiff<V1>)
                Base::getLHS().reverse(hadamard(Base::getRHS().values(), grad));
            if constexpr (ReverseDiff<V2>)
                Base::getRHS().reverse(hadamard(Base::getLHS().values(), grad));
        }
    }

    template<Vector V, Scalar U>
    [[nodiscard]] inline auto operator*(V&& v, U&& x) noexcept requires(!CUDA<V>) {
        return VectorExpr<ExprType::Mul, V&&, U&&>(std::forward<V>(v), std::forward<U>(x));
    }

    template<Scalar U, Vector V>
    [[nodiscard]] inline auto operator*(U&& x, V&& v) noexcept requires(!CUDA<V>) {
        return std::forward<V>(v) * std::forward<U>(x);
    }

    template<Vector V1, Vector V2>
    [[nodiscard]] inline auto hadamard(V1&& v1, V2&& v2) noexcept requires(!CUDA<V1> && !CUDA<V2>) {
        return VectorExpr<ExprType::Mul, V1&&, V2&&>(std::forward<V1>(v1), std::forward<V2>(v2));
    }
}

#include "MKL/Mul.h"
