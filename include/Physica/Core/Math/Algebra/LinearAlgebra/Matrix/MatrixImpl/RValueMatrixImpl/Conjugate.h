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
        using Base = RValueMatrix<Conjugate<M>>;
    protected:
        using typename Base::T;
    private:
        const M& matrix;
    public:
        Conjugate(const M& matrix_) : matrix(matrix_) {}
        /* Getters */
        [[nodiscard]] T calc(size_t row, size_t col) const { return matrix.calc(row, col).conjugate(); }
        [[nodiscard]] size_t getRow() const noexcept { return matrix.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return matrix.getCol(); }
    };

    template<Vector V>
    class ConjugateVector<V> : public RValueVector<ConjugateVector<V>> {
        using Base = RValueVector<ConjugateVector<V>>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        const V& vec;
    public:
        explicit ConjugateVector(const V& vec_) : vec(vec_) {}
        /* Operations */
        template<Vector V1, ExecutePolicy P = Sequential>
        void assign(V1& target) const;
        /* Getters */
        [[nodiscard]] T calc(size_t index) const { return vec.calc(index).conjugate(); }
        [[nodiscard]] Tv calc_value(size_t index) const { return vec.calc_value(index).conjugate(); }
        [[nodiscard]] size_t getLength() const noexcept { return vec.getLength(); }
    };

    template<Vector V>
    template<Vector V1, ExecutePolicy P>
    void ConjugateVector<V>::assign(V1& target) const {
        for (size_t i = 0; i < vec.getLength(); ++i)
            target[i] = calc(i);
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<Conjugate<M>> : public Traits<M> {};

    template<Vector V>
    class Traits<ConjugateVector<V>> : public Traits<V> {};
}
