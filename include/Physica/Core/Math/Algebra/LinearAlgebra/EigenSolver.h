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
        using EigenvalueVector = Vector<ComplexScalar<ScalarType>, MatrixType::RowAtCompile, MatrixType::MaxRowAtCompile>;
        using EigenvectorMatrix = DenseMatrix<ComplexScalar<ScalarType>,
                                              DenseMatrixOption::Column | DenseMatrixOption::Vector,
                                              MatrixType::RowAtCompile,
                                              MatrixType::RowAtCompile>;
        const MatrixType& source;
        EigenvalueVector eigenvalues;
    public:
        EigenSolver(const LValueMatrix<MatrixType>& source_, bool computeEigenvectors);
        /* Operations */
        [[nodiscard]] const EigenvalueVector& getEigenvalues() const noexcept { return eigenvalues; }
        [[nodiscard]] const EigenvectorMatrix& getEigenvectors() const;
    };

    template<class MatrixType>
    EigenSolver<MatrixType>::EigenSolver(const LValueMatrix<MatrixType>& source_, bool computeEigenvectors)
            : source(source_.getDerived())
            , eigenvalues(source_.getRow()) {
        assert(source.getRow() == source.getColumn());
        const auto schur = RealSchur(source, computeEigenvectors);

        const auto& matrixT = schur.getMatrixT();
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
                eigenvalues[i] = ComplexScalar<ScalarType>(real, -z);
                i += 2;
            }
        }
    }
}