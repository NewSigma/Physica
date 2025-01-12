/*
 * Copyright 2022-2024 Weibo He.
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

#include "RValueMatrix.h"

namespace Physica::Core {
    template<Matrix T>
    class Conjugate<T> : public RValueMatrix<Conjugate<T>> {
        using Base = RValueMatrix<Conjugate<T>>;
    public:
        using typename Base::ScalarType;
    private:
        const T& matrix;
    public:
        Conjugate(const T& matrix_) : matrix(matrix_) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return matrix.calc(row, col).conjugate(); }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return matrix.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return matrix.getCol(); }
    };

    template<Vector T>
    class ConjugateVector<T> : public RValueVector<ConjugateVector<T>> {
        using Base = RValueVector<ConjugateVector<T>>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    private:
        const T& vec;
    public:
        explicit ConjugateVector(const T& vec_) : vec(vec_) {}
        /* Operations */
        template<Vector V, class Executor = SequentialExecutor>
        void assign(V& target) const;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const { return vec.calc(index).conjugate(); }
        [[nodiscard]] ValueType calc_value(size_t index) const { return vec.calc_value(index).conjugate(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return vec.getLength(); }
    };

    template<Vector T>
    template<Vector V, class Executor>
    void ConjugateVector<T>::assign(V& target) const {
        for (size_t i = 0; i < vec.getLength(); ++i)
            target[i] = calc(i);
    }
}

namespace Physica {
    template<Matrix T>
    class Traits<Conjugate<T>> : public Traits<T> {};

    template<Vector T>
    class Traits<ConjugateVector<T>> : public Traits<T> {};
}
