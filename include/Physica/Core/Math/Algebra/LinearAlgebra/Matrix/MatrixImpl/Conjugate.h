/*
 * Copyright 2022-2024 WeiBo He.
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

namespace Physica::Core {
    template<class MatrixType> class Conjugate;

    template<class VectorType> class ConjugateVector;

    namespace Internal {
        template<class T> class Traits;

        template<class MatrixType>
        class Traits<Conjugate<MatrixType>> : public Traits<MatrixType> {};

        template<class VectorType>
        class Traits<ConjugateVector<VectorType>> : public Traits<VectorType> {};
    }

    template<class MatrixType>
    class Conjugate : public RValueMatrix<Conjugate<MatrixType>> {
    public:
        using Base = RValueMatrix<Conjugate<MatrixType>>;
        using typename Base::ScalarType;
    private:
        const MatrixType& matrix;
    public:
        Conjugate(const RValueMatrix<MatrixType>& matrix_) : matrix(matrix_.getDerived()) {}
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return matrix.calc(row, col).conjugate(); }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return matrix.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const noexcept { return matrix.getColumn(); }
    };

    template<class VectorType>
    class ConjugateVector : public RValueVector<ConjugateVector<VectorType>> {
    public:
        using Base = RValueVector<ConjugateVector<VectorType>>;
        using typename Base::ScalarType;
    private:
        const VectorType& vec;
    public:
        explicit ConjugateVector(const RValueVector<VectorType>& vec_) : vec(vec_.getDerived()) {}
        /* Operations */
        template<class OtherVector, class Executor = SequentialExecutor>
        void assignTo(LValueVector<OtherVector>& target) const;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const { return vec.calc(index).conjugate(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return vec.getLength(); }
    };

    template<class VectorType>
    template<class OtherVector, class Executor>
    void ConjugateVector<VectorType>::assignTo(LValueVector<OtherVector>& target) const {
        for (size_t i = 0; i < vec.getLength(); ++i)
            target[i] = calc(i);
    }
}
