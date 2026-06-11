/*
 * Copyright 2026 Weibo He.
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

#include "MatrixTrig.h"

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof<M, Inverse> && instanceof_tx<typename Traits<M>::ExprType, MatrixTrig>)
    class GEMV<M, V> : public RValueVector<GEMV<M, V>> {
        using This = GEMV<M, V>;
        using Base = RValueVector<This>;
    public:
        using Base::isComplex;
    protected:
        using typename Base::T;
        using typename Base::Tc;
    private:
        decay_rvalue_t<M> inv;
        decay_rvalue_t<V> rhs;
    public:
        GEMV(M&& inv, V&& rhs);
        GEMV(const This&) = default;
        GEMV(This&&) noexcept = default;
        ~GEMV() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto& target) const;
        /* Getters */
        [[nodiscard]] size_t getLength() const { return rhs.getLength(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
    };

    template<Matrix M, Vector V> requires(instanceof<M, Inverse> && instanceof_tx<typename Traits<M>::ExprType, MatrixTrig>)
    GEMV<M, V>::GEMV(M&& inv, V&& rhs) : inv(std::forward<M>(inv)), rhs(std::forward<V>(rhs)) {}

    template<Matrix M, Vector V> requires(instanceof<M, Inverse> && instanceof_tx<typename Traits<M>::ExprType, MatrixTrig>)
    template<ExecutePolicy P>
    void GEMV<M, V>::assign(Vector auto& target) const {
        using Expr = std::remove_cvref<M>::type;
        constexpr bool Unit = Traits<Expr>::Unit;
        rhs.assign(target);

        const auto& mat = inv.getExpr();
        size_t length = getLength();
        if constexpr (Traits<Expr>::Upper) {
            if constexpr (Unit) {
                for (size_t i = length - 1; i > 0; --i)
                    target.head(i) -= mat.col(i).head(i) * target[i];
            }
            else {
                for (size_t i = length - 1; i > 0; --i) {
                    T factor = target[i] / mat.calc(i, i);
                    target[i] = factor;
                    target.head(i) -= mat.col(i).head(i) * factor;
                }
                target[0] /= mat.calc(0, 0);
            }
        }
        else {
            if constexpr (Unit) {
                for (size_t i = 0; i < length - 1; ++i)
                    target.tail(i + 1) -= mat.col(i).tail(i + 1) * target[i];
            }
            else {
                size_t i = 0;
                for (; i < length - 1; ++i) {
                    T factor = target[i] / mat.calc(i, i);
                    target[i] = factor;
                    target.tail(i + 1) -= mat.col(i).tail(i + 1) * factor;
                }
                target[i] /= mat.calc(i, i);
            }
        }
    }

    template<Matrix M, Vector V> requires(instanceof<M, Inverse> && instanceof_tx<typename Traits<M>::ExprType, MatrixTrig>)
    auto&& GEMV<M, V>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.inv);
    }

    template<Matrix M, Vector V> requires(instanceof<M, Inverse> && instanceof_tx<typename Traits<M>::ExprType, MatrixTrig>)
    auto&& GEMV<M, V>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.rhs);
    }

    template<Matrix M, Vector V> requires(instanceof<M, Inverse> && instanceof_tx<typename Traits<M>::ExprType, MatrixTrig>)
    __host__ __device__ consteval size_t GEMV<M, V>::getSizeAtCompile() noexcept {
        return std::max(std::remove_cvref_t<M>::getRowAtCompile(), std::remove_cvref_t<V>::getSizeAtCompile());
    }
}
