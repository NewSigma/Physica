/*
 * Copyright 2021 WeiBo He.
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
        const MatrixType& matrix;
    public:
        using Base = RValueMatrix<Conjugate<MatrixType>>;
        using typename Base::ScalarType;
    public:
        Conjugate(const RValueMatrix<MatrixType>& matrix_) : matrix(matrix_.getDerived()) {}
        template<class OtherMatrix>
        void assignTo(LValueMatrix<OtherMatrix>& target) const;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return matrix.calc(row, col).conjugate(); }
        [[nodiscard]] size_t getRow() const noexcept { return matrix.getRow(); }
        [[nodiscard]] size_t getColumn() const noexcept { return matrix.getColumn(); }
    };

    template<class MatrixType>
    template<class OtherMatrix>
    void Conjugate<MatrixType>::assignTo(LValueMatrix<OtherMatrix>& target) const {
        using TargetType = LValueMatrix<OtherMatrix>;
        const size_t max_i = target.getMaxMajor();
        const size_t mat_j = target.getMaxMinor();
        for (size_t i = 0; i < max_i; ++i) {
            for (size_t j = 0; j < mat_j; ++j) {
                target.getElementFromMajorMinor(i, j) = calc(TargetType::rowFromMajorMinor(i, j),
                                                             TargetType::columnFromMajorMinor(i, j));
            }
        }
    }

    template<class VectorType>
    class ConjugateVector : public RValueVector<ConjugateVector<VectorType>> {
        const VectorType& vec;
    public:
        using Base = RValueVector<ConjugateVector<VectorType>>;
        using typename Base::ScalarType;
    public:
        ConjugateVector(const RValueVector<VectorType>& vec_) : vec(vec_.getDerived()) {}
        template<class OtherVector>
        void assignTo(LValueVector<OtherVector>& target) const;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const { return vec.calc(index).conjugate(); }
        [[nodiscard]] size_t getLength() const noexcept { return vec.getLength(); }
    };

    template<class VectorType>
    template<class OtherVector>
    void ConjugateVector<VectorType>::assignTo(LValueVector<OtherVector>& target) const {
        for (size_t i = 0; i < vec.getLength(); ++i)
            target[i] = calc(i);
    }
}
