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

#include <cassert>

namespace Physica::Core {
    template<class MatrixType> class InverseMatrix;

    namespace Internal {
        template<class MatrixType>
        class Traits<InverseMatrix<MatrixType>> : public Traits<MatrixType> {};
    }

    template<class MatrixType>
    class InverseMatrix : public RValueMatrix<InverseMatrix<MatrixType>> {
        const MatrixType& matrix;
    public:
        InverseMatrix(const LValueMatrix<MatrixType>& matrix_) : matrix(matrix_.getDerived()) {
            assert(matrix.getRow() == matrix.getColumn());
        }
        template<class OtherMatrix>
        void assignTo(LValueMatrix<OtherMatrix>& inverse) const;
        /* Getters */
        [[nodiscard]] const MatrixType& getMatrix() const noexcept { return matrix; }
        [[nodiscard]] size_t getRow() const noexcept { return matrix.getRow(); }
        [[nodiscard]] size_t getColumn() const noexcept { return matrix.getRow(); }
    };

    template<class MatrixType>
    template<class OtherMatrix>
    void InverseMatrix<MatrixType>::assignTo(LValueMatrix<OtherMatrix>& inverse) const {
        const size_t order = getRow();
        const size_t order_1 = order - 1;
        MatrixType copy = matrix;
        if constexpr (DenseMatrixOption::isSameMajor<MatrixType, OtherMatrix>()) {
            inverse.getDerived().toUnitMatrix();
            for (size_t i = 0; i < order_1; ++i) {
                size_t k = i;
                while(copy.majorGet(k, i).isZero()) {
                    ++k;
                    assert(k < order);
                }
                if (k != i) {
                    copy.majorSwap(k, i);
                    inverse.majorSwap(k, i);
                }
                
                for (size_t j = i + 1; j < order; ++j) {
                    auto factor = copy.majorGet(j, i) / copy.majorGet(i, i);
                    copy.majorReduce(j, i, factor);
                    inverse.majorReduce(j, i, factor);
                }
            }

            for (size_t i = order_1; i > 0; --i) {
                size_t k = i;
                while(copy.majorGet(k, i).isZero()) {
                    --k;
                    assert(k < order);
                }
                if (k != i) {
                    copy.majorSwap(k, i);
                    inverse.majorSwap(k, i);
                }
                
                for (size_t j = 0; j < i; ++j) {
                    auto factor = copy.majorGet(j, i) / copy.majorGet(i, i);
                    copy.majorReduce(j, i, factor);
                    inverse.majorReduce(j, i, factor);
                }
            }
            for (size_t i = 0; i < order; ++i)
                inverse.majorMulScalar(i, reciprocal(copy(i, i)));
        }
        else {
            auto temp = MatrixType::unitMatrix(order);
            for (size_t i = 0; i < order_1; ++i) {
                size_t k = i;
                while(copy.majorGet(k, i).isZero()) {
                    ++k;
                    assert(k < order);
                }
                if (k != i) {
                    copy.majorSwap(k, i);
                    temp.majorSwap(k, i);
                }

                for (size_t j = i + 1; j < order; ++j) {
                    auto factor = copy.majorGet(j, i) / copy.majorGet(i, i);
                    copy.majorReduce(j, i, factor);
                    temp.majorReduce(j, i, factor);
                }
            }

            for (size_t i = order_1; i > 0; --i) {
                size_t k = i;
                while(copy.majorGet(k, i).isZero()) {
                    --k;
                    assert(k < order);
                }
                if (k != i) {
                    copy.majorSwap(k, i);
                    inverse.majorSwap(k, i);
                }

                for (size_t j = 0; j < i; ++j) {
                    auto factor = copy.majorGet(j, i) / copy.majorGet(i, i);
                    copy.majorReduce(j, i, factor);
                    temp.majorReduce(j, i, factor);
                }
            }
            for (size_t i = 0; i < order; ++i)
                temp.majorMulScalar(i, reciprocal(copy(i, i)));
            inverse = temp;
        }
    }
}