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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/RValueMatrix.h"

namespace Physica {
    template<Matrix M, Vector V> class MatPowVecProd;

    template<Matrix M>
    class MatrixPow : public RValueMatrix<MatrixPow<M>> {
        static_assert(std::is_reference_v<M>);

        using This = MatrixPow<M>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
    private:
        LazyDestroy<M> m;
        int power;
    public:
        MatrixPow(M m_, int power_);
        MatrixPow(const This&) = default;
        MatrixPow(This&&) = default;
        ~MatrixPow() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;

        template<Vector V>
        [[nodiscard]] auto operator*(this auto&&, V&& v) noexcept;
        /* Operations */
        [[nodiscard]] T calc(size_t, size_t) const { noImpl("MatrixPow::calc() is low performance and should be avoided"); }

        [[nodiscard]] auto transpose() const noexcept;
        [[nodiscard]] auto hermite() const noexcept;
        /* Getters */
        [[nodiscard]] const auto& getMatrix() const noexcept { return m; }
        [[nodiscard]] size_t getRow() const noexcept { return m.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return m.getCol(); }
        [[nodiscard]] int getPower() const noexcept { return power; }
    };

    template<Matrix M>
    MatrixPow<M>::MatrixPow(M m_, int power_) : m(std::forward<M>(m_)), power(power_) {}

    template<Matrix M>
    auto MatrixPow<M>::transpose() const noexcept {
        return pow(m.transpose(), power);
    }

    template<Matrix M>
    auto MatrixPow<M>::hermite() const noexcept {
        return pow(m.hermite(), power);
    }

    template<Matrix M>
    template<Vector V>
    auto MatrixPow<M>::operator*(this auto&& self, V&& v) noexcept {
        using Self = decltype(self);
        return MatPowVecProd<Self, V&&>(std::forward<Self>(self), std::forward<V>(v));
    }

    template<Matrix M>
    [[nodiscard]] auto pow(M&& m, int power) noexcept {
        return MatrixPow<M&&>(std::forward<M>(m), power);
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<MatrixPow<M>> : public Traits<std::remove_cvref_t<M>> {
    public:
        constexpr static int Option = MatrixOption::AnyMajor;
    };
}

#include "MatPowVecProd.h"
