/*
 * Copyright 2025-2026 Weibo He.
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

#include "Inverse.h"

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof_tx<Inverse, M> && instanceof_tx<DenseLU, typename Traits<M>::ExprType>)
    class GEMV<M, V> : public RValueVector<GEMV<M, V>> {
        using This = GEMV<M, V>;
        using Base = RValueVector<This>;
        using typename Base::T;

        LazyDestroy<M> m;
        LazyDestroy<V> v;
    public:
        GEMV(M&& m, V&& v);
        GEMV(const This&) = default;
        GEMV(This&&) noexcept = default;
        ~GEMV() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto& target) const;
        void assign_base(Vector auto& target) const;
        void assign_mkl(Vector auto& target) const;

        [[nodiscard]] T calc(size_t) const { noImpl(__func__); }

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return m.getRow(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
    };

    template<Matrix M, Vector V> requires(instanceof_tx<Inverse, M> && instanceof_tx<DenseLU, typename Traits<M>::ExprType>)
    GEMV<M, V>::GEMV(M&& m, V&& v) : m(std::forward<M>(m)), v(std::forward<V>(v)) {}

    template<Matrix M, Vector V> requires(instanceof_tx<Inverse, M> && instanceof_tx<DenseLU, typename Traits<M>::ExprType>)
    template<ExecutePolicy P>
    void GEMV<M, V>::assign(Vector auto& target) const {
        if constexpr (Internal::EnableMKL<V, decltype(target)>::value)
            assign_mkl(target);
        else
            assign_base(target);
    }

    template<Matrix M, Vector V> requires(instanceof_tx<Inverse, M> && instanceof_tx<DenseLU, typename Traits<M>::ExprType>)
    void GEMV<M, V>::assign_base(Vector auto& target) const {
        const auto& lu = getLHS().getDenseLU();
        if constexpr (lu.isPivot()) {
            VectorND<T> temp1 = lu.getPerm().inv() * v;
            VectorND<T> temp2 = lu.getMatrixL().inv() * temp1;
            (lu.getMatrixU().inv() * temp2).assign(target);
        }
        else {
            VectorND<T> temp = lu.getMatrixL().inv() * v;
            (lu.getMatrixU().inv() * temp).assign(target);
        }
    }

    template<Matrix M, Vector V> requires(instanceof_tx<Inverse, M> && instanceof_tx<DenseLU, typename Traits<M>::ExprType>)
    auto GEMV<M, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() + std::forward<Self>(self).getRHS().values();
    }

    template<Matrix M, Vector V> requires(instanceof_tx<Inverse, M> && instanceof_tx<DenseLU, typename Traits<M>::ExprType>)
    auto&& GEMV<M, V>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.m);
    }

    template<Matrix M, Vector V> requires(instanceof_tx<Inverse, M> && instanceof_tx<DenseLU, typename Traits<M>::ExprType>)
    auto&& GEMV<M, V>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.v);
    }
}

#ifdef PHYSICA_MKL
    #include "InverseGEMV_MKL.h"
#endif
