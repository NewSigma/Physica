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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/MatrixExpr.h"

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof_x<M, MatrixExpr>)
    class GEMV<M, V> : public RValueVector<GEMV<M, V>> {
        using This = GEMV<M, V>;
        using Base = RValueVector<This>;
        using M1 = std::remove_cvref_t<M>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        decay_rvalue_t<M> expr;
        decay_rvalue_t<V> vec;
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

        [[nodiscard]] auto values(this auto&& self) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return expr.getRow(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isFastAssign() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return std::remove_cvref_t<M>::getRowAtCompile(); }
    };

    template<Matrix M, Vector V> requires(instanceof_x<M, MatrixExpr>)
    GEMV<M, V>::GEMV(M expr_, V vec_) : expr(std::forward<M>(expr_)), vec(std::forward<V>(vec_)) {
        assert(expr.getCol() == vec.getLength());
    }

    template<Matrix M, Vector V> requires(instanceof_x<M, MatrixExpr>)
    auto GEMV<M, V>::operator-() const noexcept {
        if constexpr (M1::isFastAssign())
            return (-getLHS()) * getRHS();
        else
            return getLHS() * (-getRHS());
    }

    template<Matrix M, Vector V> requires(instanceof_x<M, MatrixExpr>)
    template<ExecutePolicy P>
    void GEMV<M, V>::assign(Vector auto& target) const {
        if constexpr (isFastAssign()) {
            constexpr ExprID ID = std::remove_cvref_t<M>::getExprID();
            if constexpr (ID == ExprID::Minus) {
                const auto& lhs = getLHS();
                const auto& rhs = getRHS();
                if constexpr (M1::isFastAssign())
                    ((-lhs) * rhs).template assign<P>(target);
                else
                    (lhs * (-rhs)).template assign<P>(target);
            }
            else if constexpr (ID == ExprID::Add)
                (expr.getLHS() * vec + expr.getRHS() * vec).template assign<P>(target);
            else if constexpr (ID == ExprID::Sub)
                (expr.getLHS() * vec - expr.getRHS() * vec).template assign<P>(target);
            else {
                static_assert(ID == ExprID::Mul, "[Error]: assign is not implemented");
                ((expr.getLHS() * vec) * expr.getRHS()).template assign<P>(target);
            }
        }
        else
            Base::assign(target);
    }

    template<Matrix M, Vector V> requires(instanceof_x<M, MatrixExpr>)
    template<ExecutePolicy P>
    void GEMV<M, V>::assign_add(Vector auto& target) const {
        if constexpr (isFastAssign()) {
            constexpr ExprID ID = expr.getExprID();
            if constexpr (ID == ExprID::Minus) {
                const auto& lhs = getLHS();
                const auto& rhs = getRHS();
                if constexpr (M1::isFastAssign())
                    ((-lhs) * rhs).template assign_add<P>(target);
                else
                    (lhs * (-rhs)).template assign_add<P>(target);
            }
            else if constexpr (ID == ExprID::Add) {
                auto expr1 = expr.getLHS() * vec + expr.getRHS() * vec;
                expr1.template assign_add<P>(target);
            }
            else if constexpr (ID == ExprID::Sub) {
                auto expr1 = expr.getLHS() * vec - expr.getRHS() * vec;
                expr1.template assign_add<P>(target);
            }
            else {
                static_assert(ID == ExprID::Mul, "[Error]: assign is not implemented");
                auto expr1 = (expr.getLHS() * vec) * expr.getRHS();
                expr1.template assign_add<P>(target);
            }
        }
        else
            Base::assign_add(target);
    }

    template<Matrix M, Vector V> requires(instanceof_x<M, MatrixExpr>)
    auto GEMV<M, V>::calc(size_t index) const -> T {
        return expr.row(index) * vec;
    }

    template<Matrix M, Vector V> requires(instanceof_x<M, MatrixExpr>)
    auto GEMV<M, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() * std::forward<Self>(self).getRHS().values();
    }

    template<Matrix M, Vector V> requires(instanceof_x<M, MatrixExpr>)
    auto&& GEMV<M, V>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.expr);
    }

    template<Matrix M, Vector V> requires(instanceof_x<M, MatrixExpr>)
    auto&& GEMV<M, V>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.vec);
    }

    template<Matrix M, Vector V> requires(instanceof_x<M, MatrixExpr>)
    __host__ __device__ consteval bool GEMV<M, V>::isFastAssign() noexcept {
        using U = Traits<M1>::RHS;
        constexpr bool isScalarU = Scalar<U>;
        constexpr ExprID ID = M1::getExprID();

        if constexpr (ID == ExprID::Minus)
            return M1::isFastAssign();
        else if constexpr (ID == ExprID::Add || ID == ExprID::Sub) {
            if constexpr (!isScalarU) {
                using M2 = Traits<M1>::LHS;
                using Add = decltype(std::declval<M2>() * std::declval<V>() + std::declval<U>() * std::declval<V>());
                return Add::isFastAssign();
            }
            return false;
        }
        else if constexpr (ID == ExprID::Mul)
            return isScalarU;
        else
            return false;
    }
}

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof_x<M, MatrixExpr>)
    class Traits<GEMV<M, V>> {
        using M1 = std::remove_cvref_t<M>;
        using V1 = std::remove_cvref_t<V>;
        static_assert(M1::getColAtCompile() == V1::getSizeAtCompile() || M1::getColAtCompile() == Dynamic || V1::getSizeAtCompile() == Dynamic,
                      "Row and col do not match in matrix product");
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename M1::ScalarType, typename V1::ScalarType>::Type;
    };
}
