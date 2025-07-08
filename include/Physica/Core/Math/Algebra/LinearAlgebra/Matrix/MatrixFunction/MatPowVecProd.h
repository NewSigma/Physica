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

#include "MatrixPow.h"

namespace Physica {
    template<Matrix M, Vector V>
    class MatPowVecProd : public RValueVector<MatPowVecProd<M, V>> {
        using This = MatPowVecProd<M, V>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    protected:
        using typename Base::Tv;
    private:
        const LazyDestroy<M> mpow;
        const LazyDestroy<V> v;
    public:
        MatPowVecProd(M&& mpow_, V&& v_);
        MatPowVecProd(const This&) = default;
        MatPowVecProd(This&&) noexcept = default;
        ~MatPowVecProd() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto& target) const;

        [[nodiscard]] ScalarType calc(size_t) const { noImpl(__func__); }
        [[nodiscard]] Tv calc_value(size_t) const { noImpl(__func__); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
        [[nodiscard]] const MatrixPow<M>& getLHS() const noexcept { return mpow; }
        [[nodiscard]] const V& getRHS() const noexcept { return v; }
    };

    template<Matrix M, Vector V>
    MatPowVecProd<M, V>::MatPowVecProd(M&& mpow_, V&& v_) : mpow(std::forward<M>(mpow_)), v(std::forward<V>(v_)) {
        assert(mpow.getCol() == v.getLength());
    }

    template<Matrix M, Vector V>
    template<ExecutePolicy P>
    void MatPowVecProd<M, V>::assign(Vector auto& target) const {
        const int power = mpow.getPower();
        if (power == 0) {
            target = v;
            return;
        }

        using V1 = std::remove_cvref_t<decltype(target)>;
        V1 buffer = mpow.getMatrix() * v;
        for (int i = 1; i < power; ++i) {
            buffer.swap(target);
            buffer = mpow.getMatrix() * target;
        }
        buffer.swap(target);
    }
}

namespace Physica {
    template<Matrix M, Vector V>
    class Traits<MatPowVecProd<M, V>> : public Traits<GEMV<M, V>> {};
}
