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

#include "GEMV.h"

namespace Physica {
    template<ExprType, class, class U> class MatrixExpr;

    template<Matrix M, Vector V> requires(instanceof_xt<MatrixExpr, M>)
    class GEMV<M, V> : public RValueVector<GEMV<M, V>> {
        using This = GEMV<M, V>;
        using Base = RValueVector<This>;
        using M1 = std::remove_cvref_t<M>;
        constexpr static ExprType Type = Traits<M1>::Type;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        LazyDestroy<M> expr;
        LazyDestroy<V> vec;
    public:
        GEMV(M expr_, V vec_);
        GEMV(const This&) = default;
        GEMV(This&&) noexcept = default;
        ~GEMV() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        [[nodiscard]] auto operator-() const noexcept;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto& target) const;
        template<ExecutePolicy P = Sequential>
        void assign_add(Vector auto& target) const;

        [[nodiscard]] T calc(size_t index) const;
        [[nodiscard]] Tv calc_value(size_t index) const;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return expr.getRow(); }
        [[nodiscard]] const auto& getLHS() const noexcept { return expr; }
        [[nodiscard]] const auto& getRHS() const noexcept { return vec; }
    };

    template<Matrix M, Vector V> requires(instanceof_xt<MatrixExpr, M>)
    GEMV<M, V>::GEMV(M expr_, V vec_) : expr(std::forward<M>(expr_)), vec(std::forward<V>(vec_)) {
        assert(expr.getCol() == vec.getLength());
    }

    template<Matrix M, Vector V> requires(instanceof_xt<MatrixExpr, M>)
    auto GEMV<M, V>::operator-() const noexcept {
        if constexpr (Traits<M1>::FastAssign)
            return (-getLHS()) * getRHS();
        else
            return getLHS() * (-getRHS());
    }

    template<Matrix M, Vector V> requires(instanceof_xt<MatrixExpr, M>)
    template<ExecutePolicy P>
    void GEMV<M, V>::assign(Vector auto& target) const {
        constexpr bool FastAssign = Traits<This>::FastAssign;
        if constexpr (FastAssign) {
            if constexpr (Type == ExprType::Minus) {
                const auto& lhs = getLHS();
                const auto& rhs = getRHS();
                if constexpr (Traits<std::remove_cvref_t<M>>::FastAssign)
                    ((-lhs) * rhs).template assign<P>(target);
                else
                    (lhs * (-rhs)).template assign<P>(target);
            }
            else if constexpr (Type == ExprType::Add)
                (expr.getLHS() * vec + expr.getRHS() * vec).template assign<P>(target);
            else if constexpr (Type == ExprType::Sub)
                (expr.getLHS() * vec - expr.getRHS() * vec).template assign<P>(target);
            else if constexpr (Type == ExprType::Mul)
                ((expr.getLHS() * vec) * expr.getRHS()).template assign<P>(target);
            else
                static_assert(!FastAssign, "[Error]: assign is not implemented");
        }
        else
            Base::assign(target);
    }

    template<Matrix M, Vector V> requires(instanceof_xt<MatrixExpr, M>)
    template<ExecutePolicy P>
    void GEMV<M, V>::assign_add(Vector auto& target) const {
        constexpr bool FastAssign = Traits<This>::FastAssign;
        if constexpr (FastAssign) {
            if constexpr (Type == ExprType::Minus) {
                const auto& lhs = getLHS();
                const auto& rhs = getRHS();
                if constexpr (Traits<std::remove_cvref_t<M>>::FastAssign)
                    ((-lhs) * rhs).template assign_add<P>(target);
                else
                    (lhs * (-rhs)).template assign_add<P>(target);
            }
            else if constexpr (Type == ExprType::Add) {
                auto expr1 = expr.getLHS() * vec + expr.getRHS() * vec;
                expr1.template assign_add<P>(target);
            }
            else if constexpr (Type == ExprType::Sub) {
                auto expr1 = expr.getLHS() * vec - expr.getRHS() * vec;
                expr1.template assign_add<P>(target);
            }
            else if constexpr (Type == ExprType::Mul) {
                auto expr1 = (expr.getLHS() * vec) * expr.getRHS();
                expr1.template assign_add<P>(target);
            }
            else
                static_assert(!FastAssign, "[Error]: assign is not implemented");
        }
        else
            Base::assign_add(target);
    }

    template<Matrix M, Vector V> requires(instanceof_xt<MatrixExpr, M>)
    auto GEMV<M, V>::calc(size_t index) const -> T {
        return expr.row(index) * vec;
    }

    template<Matrix M, Vector V> requires(instanceof_xt<MatrixExpr, M>)
    auto GEMV<M, V>::calc_value(size_t index) const -> Tv {
        return expr.row(index).values() * vec.values();
    }
}

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof_xt<MatrixExpr, M>)
    class Traits<GEMV<M, V>> {
        using M1 = std::remove_cvref_t<M>;
        using V1 = std::remove_cvref_t<V>;
        static_assert(M1::ColAtCompile == V1::SizeAtCompile || M1::ColAtCompile == Dynamic || V1::SizeAtCompile == Dynamic,
                      "Row and col do not match in matrix product");

        constexpr static bool calcFastAssign() {
            using U = Traits<M1>::RHS;
            constexpr bool isScalarU = Scalar<U>;
            constexpr ExprType Type = Traits<M1>::Type;

            if constexpr (Type == ExprType::Minus)
                return Traits<M1>::FastAssign;
            else if constexpr (Type == ExprType::Add || Type == ExprType::Sub) {
                if constexpr (!isScalarU) {
                    using M2 = Traits<M1>::LHS;
                    using Add = decltype(std::declval<M2>() * std::declval<V>() + std::declval<U>() * std::declval<V>());
                    return Traits<Add>::FastAssign;
                }
                return false;
            }
            else if constexpr (Type == ExprType::Mul)
                return isScalarU;
            else
                return false;
        }
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename M1::ScalarType, typename V1::ScalarType>::Type;
        constexpr static size_t SizeAtCompile = M1::RowAtCompile;

        constexpr static bool FastAssign = calcFastAssign();
        constexpr static bool FastPacket = false;
    };
}
