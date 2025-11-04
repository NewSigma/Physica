/*
 * Copyright 2022-2025 Weibo He.
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
        static_assert(M::isComplex, "[Error]: Unnecessary conjugate call on real matrix");

        using This = Conjugate<M>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
    private:
        const M& matrix;
    public:
        Conjugate(const M& matrix_) : matrix(matrix_) {}
        Conjugate(const This&) = delete;
        Conjugate(This&&) = delete;
        ~Conjugate() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const { return matrix.calc(row, col).conjugate(); }

        [[nodiscard]] const M& conjugate() const noexcept { return matrix; }
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return matrix.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return matrix.getCol(); }
    };
}

namespace Physica {
    template<Matrix M>
    class Traits<Conjugate<M>> : public Traits<M> {};
}
