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

#include "../MatrixExpr.h"

namespace Physica {
    template<Matrix M, Vector V>
    class MatExprVecProd : public RValueVector<MatExprVecProd<M, V>> {
        using This = MatExprVecProd<M, V>;
        using Base = RValueVector<This>;
        using M1 = std::remove_cvref_t<M>;
        constexpr static ExprType Type = Traits<M1>::Type;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        const LazyDestroy<M> expr;
        const LazyDestroy<V> vec;
    public:
        MatExprVecProd(M&& expr_, V&& vec_);
        MatExprVecProd(const This&) = default;
        MatExprVecProd(This&&) noexcept = default;
        ~MatExprVecProd() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        [[nodiscard]] auto operator-() const noexcept;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        inline void assign(Vector auto& target) const;
        template<ExecutePolicy P = Sequential>
        inline void assign_add(Vector auto& target) const;

        [[nodiscard]] inline T calc(size_t index) const;
        [[nodiscard]] inline Tv calc_value(size_t index) const;
        /* Getters */
        [[nodiscard]] size_t getLength() const { return expr.getRow(); }
        [[nodiscard]] const auto& getLHS() const noexcept { return expr; }
        [[nodiscard]] const auto& getRHS() const noexcept { return vec; }
    };

    template<Matrix M, Vector V>
    MatExprVecProd<M, V>::MatExprVecProd(M&& expr_, V&& vec_) : expr(std::forward<M>(expr_)), vec(std::forward<V>(vec_)) {
        assert(expr.getCol() == vec.getLength());
    }

    template<Matrix M, Vector V>
    auto MatExprVecProd<M, V>::operator-() const noexcept {
        if constexpr (Traits<M1>::FastAssign)
            return (-getLHS()) * getRHS();
        else
            return getLHS() * (-getRHS());
    }

    template<Matrix M, Vector V>
    template<ExecutePolicy P>
    inline void MatExprVecProd<M, V>::assign(Vector auto& target) const {
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
                target.template operator=<P>(expr.getLHS() * vec + expr.getRHS() * vec);
            else if constexpr (Type == ExprType::Sub)
                target.template operator=<P>(expr.getLHS() * vec - expr.getRHS() * vec);
            else if constexpr (Type == ExprType::Mul)
                target.template operator=<P>((expr.getLHS() * vec) * expr.getRHS());
            else
                static_assert(!FastAssign, "[Error]: assign is not implemented");
        }
        else
            Base::assign(target);
    }

    template<Matrix M, Vector V>
    template<ExecutePolicy P>
    inline void MatExprVecProd<M, V>::assign_add(Vector auto& target) const {
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

    template<Matrix M, Vector V>
    inline auto MatExprVecProd<M, V>::calc(size_t index) const -> T {
        return expr.row(index) * vec;
    }

    template<Matrix M, Vector V>
    inline auto MatExprVecProd<M, V>::calc_value(size_t index) const -> Tv {
        return expr.row(index).values() * vec.values();
    }
}

namespace Physica {
    template<Matrix M, Vector V>
    class Traits<MatExprVecProd<M, V>> {
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
