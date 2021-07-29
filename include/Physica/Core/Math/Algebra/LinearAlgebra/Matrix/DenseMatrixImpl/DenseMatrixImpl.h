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

namespace Physica::Core {
    template<class T, int type, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    template<Utils::ExpressionType expType, class T1, class T2>
    DenseMatrix<T, type, Row, Column, MaxRow, MaxColumn>::DenseMatrix(DenseMatrixExpression<expType, T1, T2> exp)
            : DenseMatrix(exp.getRow(), exp.getColumn()) {
        for (size_t i = 0; i < Base::getMaxMajor(); i++) {
            for (size_t j = 0; j < Base::getMaxMinor(); ++j) {
                Base::getElementFromMajorMinor(i, j) = exp(Base::rowFromMajorMinor(i, j), Base::columnFromMajorMinor(i, j));
            }
        }
    }

    template<class T, int type, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    template<class OtherMatrix>
    DenseMatrix<T, type, Row, Column, MaxRow, MaxColumn>::DenseMatrix(const DenseMatrixBase<OtherMatrix>& mat)
            : DenseMatrix(mat.getRow(), mat.getColumn()) {
        mat.assignTo(*this);
    }

    template<class T, int type, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    template<class T1, class T2>
    DenseMatrix<T, type, Row, Column, MaxRow, MaxColumn>::DenseMatrix(MatrixProduct<T1, T2> pro) : DenseMatrix(pro.getRow(), pro.getColumn()) {
        for (size_t i = 0; i < Base::getMaxMajor(); i++) {
            for (size_t j = 0; j < Base::getMaxMinor(); ++j) {
                Base::getElementFromMajorMinor(i, j) = pro(Base::rowFromMajorMinor(i, j), Base::columnFromMajorMinor(i, j));
            }
        }
    }

    template<class T, int type, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    template<class MatrixIn>
    DenseMatrix<T, type, Row, Column, MaxRow, MaxColumn>::DenseMatrix(LUDecomposition<MatrixIn> lu)
            : DenseMatrix(lu.getOrder(), lu.getOrder()) {
        const size_t rank = lu.getOrder();
        (*this) = lu.getMatrix();
        for (size_t i = 0; i < rank; ++i)
            lu.decompositionColumn((*this), i);
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
