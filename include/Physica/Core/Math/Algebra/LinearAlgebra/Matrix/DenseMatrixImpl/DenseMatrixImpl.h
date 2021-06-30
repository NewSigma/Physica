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
namespace Physica::Core {
    template<class T, int type, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    template<class MatrixIn>
    DenseMatrix<T, type, Row, Column, MaxRow, MaxColumn>::DenseMatrix(LUDecomposition<MatrixIn> lu)
            : DenseMatrix(lu.getOrder(), lu.getOrder()) {
        const size_t rank = lu.getOrder();
        (*this) = lu.getMatrix();
        for (size_t i = 0; i < rank; ++i)
            lu.decompositionColumn((*this), i);
    }
    /**
     * Implemented the square method
     */
    template<class T, int type, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    template<class MatrixIn>
    DenseMatrix<T, type, Row, Column, MaxRow, MaxColumn>::DenseMatrix(Cholesky<MatrixIn> cholesky)
            : DenseMatrix(cholesky.getOrder(), cholesky.getOrder()) {
        typedef DenseMatrix<T, type, Row, Column, MaxRow, MaxColumn> MatrixOut;
        const size_t order = cholesky.getOrder();
        const MatrixIn& matrix = cholesky.getMatrix();
        auto matrixIte = Base::begin();
        auto constMatrixIte = matrix.cbegin();

        auto elementIte = this->ebegin(matrixIte);
        auto constElementIte = matrix.cebegin(constMatrixIte);
        /* Handle first vector */ {
            const auto diag = sqrt(*constElementIte);
            *elementIte = diag;
            for (size_t i = 1; i < order; ++i) {
                ++elementIte;
                ++constElementIte;
                *elementIte = *constElementIte / diag;
            }
        }
        /* Handle other vectors */ {
            for (size_t i = 1; i < order; ++i) {
                MatrixOut::updateIterator(matrixIte, elementIte);
                MatrixIn::updateIterator(constMatrixIte, constElementIte);
                size_t j;
                for (j = 0; j < i; ++j) {
                    *elementIte = 0;
                    ++elementIte;
                    ++constElementIte;
                }

                T diag(*constElementIte);
                /* j == i */ {
                    for (size_t k = 0; k < i; ++k) {
                        if constexpr (DenseMatrixType::isColumnMatrix<MatrixOut>())
                            diag -= square((*this)(i, k));
                        else
                            diag -= square((*this)(k, i));
                    }
                    diag = sqrt(diag);
                    *elementIte = diag;
                    ++j;
                }

                for (; j < order; ++j) {
                    ++elementIte;
                    ++constElementIte;
                    T temp(*constElementIte);
                    for (size_t k = 0; k < j; ++k) {
                        if constexpr (DenseMatrixType::isColumnMatrix<MatrixOut>())
                            diag -= (*this)(i, k) * (*this)(j, k);
                        else
                            diag -= (*this)(k, i) * (*this)(k, j);
                    }
                    *elementIte = temp / diag;
                }
            }
        }
    }

    template<class T, int type, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    template<class MatrixIn>
    DenseMatrix<T, type, Row, Column, MaxRow, MaxColumn>::DenseMatrix(InverseMatrix<MatrixIn> inverse)
            : DenseMatrix(inverse.getOrder(), inverse.getOrder()) {
        using MatrixOut = DenseMatrix<T, type, Row, Column, MaxRow, MaxColumn>;
        MatrixIn matrixIn(inverse.getMatrix());
        const size_t order = inverse.getOrder();
        const size_t order_1 = order - 1;
        if constexpr (DenseMatrixType::isSameMajor<MatrixIn, MatrixOut>()) {
            Base::toUnitMatrix();
            for (size_t i = 0; i < order_1; ++i) {
                size_t k = i;
                while(matrixIn.majorGet(k, i).isZero()) {
                    ++k;
                    assert(k < order);
                }
                if (k != i) {
                    matrixIn.majorSwap(k, i);
                    Base::majorSwap(k, i);
                }
                
                for (size_t j = i + 1; j < order; ++j) {
                    auto factor = matrixIn.majorGet(j, i) / matrixIn.majorGet(i, i);
                    matrixIn.majorReduce(j, i, factor);
                    Base::majorReduce(j, i, factor);
                }
            }

            for (size_t i = order_1; i > 0; --i) {
                size_t k = i;
                while(matrixIn.majorGet(k, i).isZero()) {
                    --k;
                    assert(k < order);
                }
                if (k != i) {
                    matrixIn.majorSwap(k, i);
                    Base::majorSwap(k, i);
                }
                
                for (size_t j = 0; j < i; ++j) {
                    auto factor = matrixIn.majorGet(j, i) / matrixIn.majorGet(i, i);
                    matrixIn.majorReduce(j, i, factor);
                    Base::majorReduce(j, i, factor);
                }
            }
            for (size_t i = 0; i < order; ++i)
                Base::majorMulScalar(i, reciprocal(matrixIn(i, i)));
        }
        else {
            auto temp = MatrixIn::unitMatrix(order);
            for (size_t i = 0; i < order_1; ++i) {
                size_t k = i;
                while(matrixIn.majorGet(k, i).isZero()) {
                    ++k;
                    assert(k < order);
                }
                if (k != i) {
                    matrixIn.majorSwap(k, i);
                    temp.majorSwap(k, i);
                }

                for (size_t j = i + 1; j < order; ++j) {
                    auto factor = matrixIn.majorGet(j, i) / matrixIn.majorGet(i, i);
                    matrixIn.majorReduce(j, i, factor);
                    temp.majorReduce(j, i, factor);
                }
            }

            for (size_t i = order_1; i > 0; --i) {
                size_t k = i;
                while(matrixIn.majorGet(k, i).isZero()) {
                    --k;
                    assert(k < order);
                }
                if (k != i) {
                    matrixIn.majorSwap(k, i);
                    Base::majorSwap(k, i);
                }

                for (size_t j = 0; j < i; ++j) {
                    auto factor = matrixIn.majorGet(j, i) / matrixIn.majorGet(i, i);
                    matrixIn.majorReduce(j, i, factor);
                    temp.majorReduce(j, i, factor);
                }
            }
            for (size_t i = 0; i < order; ++i)
                temp.majorMulScalar(i, reciprocal(matrixIn(i, i)));
            (*this) = temp;
        }
    }

    template<class T, int type, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    DenseMatrix<T, type, Row, Column, MaxRow, MaxColumn> DenseMatrix<T, type, Row, Column, MaxRow, MaxColumn>::unitMatrix(size_t order) {
        DenseMatrix result(order, order);
        result.toUnitMatrix();
        return result;
    }

    template<class T, int type, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    std::ostream& operator<<(std::ostream& os, const DenseMatrix<T, type, Row, Column, MaxRow, MaxColumn>& mat) {
        const size_t column = mat.getColumn();
        const size_t row = mat.getRow();
        size_t width = 0;
        /* Get max width */ {
            for (size_t c = 0; c < column; ++c) {
                for (size_t r = 0; r < row; ++ r) {
                    std::stringstream stream{};
                    stream.copyfmt(os);
                    stream << mat(r, c);
                    width = std::max(width, stream.str().length());
                }
            }
        }
        /* Output */ {
            for (size_t c = 0; c < column; ++c) {
                os.width(width);
                os << mat(0, c) << ' ';
            }
            for (size_t r = 1; r < row; ++r) {
                os << '\n';
                for (size_t c = 0; c < column; ++c) {
                    os.width(width);
                    os << mat(r, c) << ' ';
                }
            }
        }
        return os;
    }
}
