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

#include "Matrix/MatrixDecomposition/RealSchur.h"
#include "Physica/Core/MultiPrecision/ComplexScalar.h"

namespace Physica::Core {
    /**
     * References:
     * [1] Golub, GeneH. Matrix computations = 矩阵计算 / 4th edition[M]. 人民邮电出版社, 2014.
     * [2] Eigen https://eigen.tuxfamily.org/
     */
    template<class MatrixType>
    class EigenSolver {
        using ScalarType = typename MatrixType::ScalarType;
    private:
        using EigenvectorStorage = DenseMatrix<ScalarType,
                                               DenseMatrixOption::Column | DenseMatrixOption::Vector,
                                               MatrixType::RowAtCompile,
                                               MatrixType::RowAtCompile>;
    public:
        using EigenvalueVector = Vector<ComplexScalar<ScalarType>, MatrixType::RowAtCompile, MatrixType::MaxRowAtCompile>;
        using EigenvectorMatrix = DenseMatrix<ComplexScalar<ScalarType>,
                                              DenseMatrixOption::Column | DenseMatrixOption::Vector,
                                              MatrixType::RowAtCompile,
                                              MatrixType::RowAtCompile>;
    private:
        const MatrixType& source;
        EigenvalueVector eigenvalues;
        EigenvectorStorage eigenvectorStorage;
    public:
        EigenSolver(const LValueMatrix<MatrixType>& source_, bool computeEigenvectors);
        /* Operations */
        [[nodiscard]] const EigenvalueVector& getEigenvalues() const noexcept { return eigenvalues; }
        [[nodiscard]] EigenvectorMatrix getEigenvectors() const;
    };

    template<class MatrixType>
    EigenSolver<MatrixType>::EigenSolver(const LValueMatrix<MatrixType>& source_, bool computeEigenvectors)
            : source(source_.getDerived())
            , eigenvalues(source_.getRow())
            , eigenvectorStorage(source_.getRow(), source_.getRow()) {
        assert(source.getRow() == source.getColumn());
        auto schur = RealSchur(source, computeEigenvectors);

        auto& matrixT = schur.getMatrixT();
        const size_t order = source.getRow();
        for (size_t i = 0; i < order;) {
            if (i == order - 1 || matrixT(i + 1, i).isZero()) {
                eigenvalues[i] = matrixT(i, i);
                i += 1;
            }
            else {
                const ScalarType p = ScalarType(0.5) * (matrixT(i, i) - matrixT(i + 1, i + 1));
                ScalarType z;
                /* Referenced from eigen, to avoid overflow */ {
                    ScalarType t0 = matrixT(i + 1, i);
                    ScalarType t1 = matrixT(i, i + 1);
                    const ScalarType max = std::max(abs(p), std::max(abs(t0), abs(t1)));
                    const ScalarType inv_max = reciprocal(max);
                    t0 *= inv_max;
                    t1 *= inv_max;
                    ScalarType p0 = p * inv_max;
                    z = max * sqrt(abs(square(p0) + t0 * t1));
                }
                const ScalarType real = p + matrixT(i + 1, i + 1);
                eigenvalues[i] = ComplexScalar<ScalarType>(real, z);
                eigenvalues[i + 1] = ComplexScalar<ScalarType>(real, -z);
                i += 2;
            }
        }

        if (computeEigenvectors) {
            for (size_t i = order - 1; i < order; --i) {
                auto block = matrixT.topLeftCorner(i + 1, i + 1);
                if (eigenvalues[i].getImag().isZero()) {
                    auto col = block.col(i);
                    col[i] = ScalarType::One();
                    for (size_t j = i - 1; j < i; --j) {
                        if (eigenvalues[j].getImag().isZero()) {
                            auto row = block.row(j);
                            col[j] = (col.tail(j + 1) * row.tail(j + 1)) / (eigenvalues[i].getReal() - eigenvalues[j].getReal());
                        }
                        else {
                            auto row1 = block.row(j - 1);
                            auto row2 = block.row(j);
                            auto tail = col.tail(j + 1);
                            ScalarType dot1 = tail * row1.tail(j + 1);
                            ScalarType dot2 = tail * row2.tail(j + 1);
                            ScalarType inv_determinate = reciprocal(square(eigenvalues[j].getReal() - eigenvalues[i].getReal()) + square(eigenvalues[j].getImag()));
                            col[j - 1] = (dot1 * block(j, j) - dot2 * block(j - 1, j)) * inv_determinate;
                            col[j] = (dot2 * block(j - 1, j - 1) - dot1 * block(j, j - 1)) * inv_determinate;
                            --j;
                        }
                    }
                }
                else {
                    assert(eigenvalues[i].getImag().isNegative());
                    auto col1 = block.col(i - 1);
                    auto col2 = block.col(i);
                    //Referenced from eigen, to ensure numerical stable
                    if (abs(matrixT(i, i - 1)) > abs(matrixT(i - 1, i))) {
                        auto temp = reciprocal(matrixT(i, i - 1));
                        col1[i - 1] = eigenvalues[i].getImag() * temp;
                        col2[i - 1] = (eigenvalues[i].getReal() - matrixT(i, i)) * temp;
                    }
                    else {
                        ComplexScalar<ScalarType> c = ComplexScalar<ScalarType>(ScalarType::Zero(), -matrixT(i - 1, i)) /
                                                      ComplexScalar<ScalarType>(matrixT(i - 1, i - 1) - eigenvalues[i].getReal(), eigenvalues[i].getImag());
                        col1[i - 1] = c.getReal();
                        col2[i - 1] = c.getImag();
                    }
                    col1[i] = ScalarType::Zero();
                    col2[i] = ScalarType::One();

                    for (size_t j = i - 2; j < i; --j) {
                        if (eigenvalues[j].getImag().isZero()) {
                            auto row = block.row(j);
                            auto tail = row.tail(j + 1);
                            const ScalarType dot1 = -(tail * col1.tail(j + 1));
                            const ScalarType dot2 = -(tail * col2.tail(j + 1));
                            const ScalarType a = block(j, j) - eigenvalues[i].getReal();
                            const ScalarType b = eigenvalues[i].getImag();
                            const ScalarType inv_denominator = reciprocal(square(a) + square(b));
                            col1[j] = (a * dot1 + b * dot2) * inv_denominator;
                            col2[j] = (a * dot2 - b * dot1) * inv_denominator;
                        }
                        else {
                            auto row1 = block.row(j - 1);
                            auto tail1 = row1.tail(j + 1);
                            const ScalarType dot11 = tail1 * col1.tail(j + 1);
                            const ScalarType dot12 = tail1 * col2.tail(j + 1);
                            auto row2 = block.row(j);
                            auto tail2 = row2.tail(j + 1);
                            const ScalarType dot21 = tail2 * col1.tail(j + 1);
                            const ScalarType dot22 = tail2 * col2.tail(j + 1);

                            const ScalarType x = matrixT(j - 1, j);
                            const ScalarType y = matrixT(j, j - 1);
                            const ScalarType vr = square(eigenvalues[j].getReal() - eigenvalues[i].getReal()) + square(eigenvalues[j].getImag()) - square(eigenvalues[i].getImag());
                            const ScalarType vi = ScalarType::Two() * (eigenvalues[j].getReal() - eigenvalues[i].getReal()) * eigenvalues[i].getImag();
                            
                            const ScalarType temp1 = matrixT(j, j) - eigenvalues[i].getReal();

                            ComplexScalar<ScalarType> c = ComplexScalar<ScalarType>(x * dot21 - temp1 * dot11 + eigenvalues[i].getImag() * dot12,
                                                                                    x * dot22 - temp1 * dot12 - eigenvalues[i].getImag() * dot11)
                                                          / ComplexScalar(vr,vi);
                            matrixT(j - 1, i - 1) = c.getReal();
                            matrixT(j - 1, i) = c.getImag();
                            if (abs(x) > (abs(temp1) + abs(eigenvalues[i].getImag()))) {
                                const ScalarType temp2 = matrixT(j - 1, j - 1) - eigenvalues[i].getReal();
                                matrixT(j, i - 1) = (-dot11 - temp2 * matrixT(j - 1, i - 1) + eigenvalues[i].getImag() * matrixT(j - 1, i)) / x;
                                matrixT(j, i) = (-dot12 - temp2 * matrixT(j - 1, i) - eigenvalues[i].getImag() * matrixT(j - 1, i - 1)) / x;
                            }
                            else {
                                c = ComplexScalar(-dot21 - y * matrixT(j - 1, i - 1), -dot22 - y * matrixT(j - 1, i)) / ComplexScalar(temp1, eigenvalues[i].getImag());
                                matrixT(j, i - 1) = c.getReal();
                                matrixT(j, i) = c.getImag();
                            }
                            --j;
                        }
                    }
                    --i;
                }
            }
            const auto& matrixU = schur.getMatrixU();
            for (size_t i = 0; i < order; ++i) {
                auto topRows = matrixT.topRows(i + 1);
                auto toCol = eigenvectorStorage.col(i);
                toCol = (matrixU.leftCols(i + 1) * topRows.col(i));
            }
        }
    }

    template<class MatrixType>
    typename EigenSolver<MatrixType>::EigenvectorMatrix EigenSolver<MatrixType>::getEigenvectors() const {
        const size_t order = eigenvalues.getLength();
        EigenvectorMatrix result = EigenvectorMatrix(order, order);
        for (size_t i = 0; i < order; ++i) {
            if (eigenvalues[i].getImag().isZero()) {
                auto toCol = result.col(i);
                toCol.asVector() = eigenvectorStorage.col(i);
            }
            else {
                auto toCol1 = result.col(i);
                auto toCol2 = result.col(i + 1);
                auto fromCol1 = eigenvectorStorage.col(i);
                auto fromCol2 = eigenvectorStorage.col(i + 1);
                for (size_t j = 0; j < order; ++j) {
                    toCol1[j] = ComplexScalar<ScalarType>(fromCol1[j], fromCol2[j]);
                    toCol2[j] = ComplexScalar<ScalarType>(fromCol1[j], -fromCol2[j]);
                }
                ++i;
            }
        }
        return result;
    }
}
