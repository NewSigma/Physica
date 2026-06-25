/*
 * Copyright 2022-2026 Weibo He.
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

#include "../RValueMatrix.h"

namespace Physica {
    template<Matrix M>
    class Conjugate<M> : public RValueMatrix<Conjugate<M>> {
        using This = Conjugate<M>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<M> m;
    public:
        Conjugate(M&& m);
        Conjugate(const This&) = default;
        Conjugate(This&&) = default;
        ~Conjugate() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const { return m.calc(row, col).conjugate(); }

        [[nodiscard]] auto&& conjugate(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return m.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return m.getCol(); }
        [[nodiscard]] size_t getOrder() const noexcept { return m.getOrder(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSquare() noexcept;
    };

    template<Matrix M>
    Conjugate<M>::Conjugate(M&& m) : m(std::forward<M>(m)) {
        static_assert(m.isComplex(), "[Error]: Unnecessary conjugate call on real matrix");
    }

    template<Matrix M>
    auto&& Conjugate<M>::conjugate(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.m);
    }

    template<Matrix M>
    __host__ __device__ consteval bool Conjugate<M>::isStaticSquare() noexcept {
        return std::remove_cvref_t<M>::isStaticSquare();
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<Conjugate<M>> : public Traits<M> {};
}
