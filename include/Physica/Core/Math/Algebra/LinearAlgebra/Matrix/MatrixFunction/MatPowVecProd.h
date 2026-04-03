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

#include "MatrixPow.h"

namespace Physica {
    template<Matrix M, Vector V>
    class MatPowVecProd : public RValueVector<MatPowVecProd<M, V>> {
        static_assert(std::is_reference_v<M> && std::is_reference_v<V>);

        using This = MatPowVecProd<M, V>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        LazyDestroy<M> mpow;
        LazyDestroy<V> v;
    public:
        MatPowVecProd(M mpow_, V v_);
        MatPowVecProd(const This&) = default;
        MatPowVecProd(This&&) noexcept = default;
        ~MatPowVecProd() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto& target) const __restrict;

        [[nodiscard]] T calc(size_t) const { noImpl(__func__); }

        [[nodiscard]] auto values(this auto&& self) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
    };

    template<Matrix M, Vector V>
    MatPowVecProd<M, V>::MatPowVecProd(M mpow_, V v_) : mpow(std::forward<M>(mpow_)), v(std::forward<V>(v_)) {
        assert(mpow.getCol() == v.getLength());
    }

    template<Matrix M, Vector V>
    template<ExecutePolicy P>
    void MatPowVecProd<M, V>::assign(Vector auto& target) const __restrict {
        target.assert_assign(v);
        const int power = mpow.getPower();
        if (power == 0) {
            v.template assign<P>(target);
            return;
        }

        auto buffer = DenseVector<T, Base::SizeAtCompile>(getLength());
        (mpow.getMatrix() * v).template assign<P>(target);
        for (int i = 1; i < power; ++i) {
            (mpow.getMatrix() * target).template assign<P>(buffer);
            buffer.swap(target);
        }
    }

    template<Matrix M, Vector V>
    auto MatPowVecProd<M, V>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS().values() * std::forward<Self>(self).getRHS().values();
    }

    template<Matrix M, Vector V>
    auto&& MatPowVecProd<M, V>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mpow);
    }

    template<Matrix M, Vector V>
    auto&& MatPowVecProd<M, V>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.v);
    }
}

namespace Physica {
    template<Matrix M, Vector V>
    class Traits<MatPowVecProd<M, V>> : public Traits<GEMV<M, V>> {};
}
