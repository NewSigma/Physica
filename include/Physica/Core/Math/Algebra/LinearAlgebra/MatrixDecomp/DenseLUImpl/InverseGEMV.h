/*
 * Copyright 2025 Weibo He.
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
    template<Matrix M, Vector V> requires(instanceof_tx<Inverse, M> && requires { std::declval<M>().getDenseLU(); })
    class GEMV<M, V> : public RValueVector<GEMV<M, V>> {
        using This = GEMV<M, V>;
        using Base = RValueVector<This>;

        using Base::isComplex;
        using typename Base::T;
        using typename Base::Tm;

        LazyDestroy<M> m;
        LazyDestroy<V> v;
    public:
        GEMV(M m, V v);
        GEMV(const This&) = default;
        GEMV(This&&) noexcept = default;
        ~GEMV() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto& target) const;
        void assign_mkl(Vector auto& target) const;

        [[nodiscard]] T calc(size_t) const { noImpl(__func__); }

        [[nodiscard]] auto values() const noexcept { return m.values() * v.values(); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return m.getRow(); }
        [[nodiscard]] const auto& getLHS() const noexcept { return m; }
        [[nodiscard]] const auto& getRHS() const noexcept { return v; }
    };

    template<Matrix M, Vector V> requires(instanceof_tx<Inverse, M> && requires { std::declval<M>().getDenseLU(); })
    GEMV<M, V>::GEMV(M m, V v) : m(std::forward<M>(m)), v(std::forward<V>(v)) {}

    template<Matrix M, Vector V> requires(instanceof_tx<Inverse, M> && requires { std::declval<M>().getDenseLU(); })
    template<ExecutePolicy P>
    void GEMV<M, V>::assign(Vector auto& target) const {
        target.assert_assign(*this);
        if constexpr (HasMKL())
            assign_mkl(target);
        else
            noImpl(__func__);
    }
}

#ifdef PHYSICA_MKL
    #include "InverseGEMV_MKL.h"
#endif
