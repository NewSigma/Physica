/*
 * Copyright 2021-2022 WeiBo He.
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

#include "Hessenburg.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Givens.h"
#include "Physica/Core/Exception/BadConvergenceException.h"

namespace Physica::Core {
    /**
     * References:
     * [1] Golub, GeneH. Matrix computations = 矩阵计算 / 4th edition[M]. 人民邮电出版社, 2014.
     * [2] Eigen https://eigen.tuxfamily.org/
     */
    template<class MatrixType>
    class Schur {
    public:
        using ScalarType = typename MatrixType::ScalarType;
        using RealType = typename ScalarType::RealType;
        using ComplexType = ComplexScalar<RealType>;
        using WorkingMatrix = typename MatrixType::ColMatrix;
        constexpr static size_t matItePerCol = 40; //Reference to Eigen
    private:
        WorkingMatrix matrixT;
        WorkingMatrix matrixU;
        bool computeMatrixU;
    public:
        template<class OtherMatrix>
        Schur(const RValueMatrix<OtherMatrix>& source, bool computeMatrixU_ = false);
        /* Getters */
        [[nodiscard]] WorkingMatrix& getMatrixT() noexcept { return matrixT; }
        [[nodiscard]] WorkingMatrix& getMatrixU() noexcept { assert(computeMatrixU); return matrixU; }
        [[nodiscard]] const WorkingMatrix& getMatrixT() const noexcept { return matrixT; }
        [[nodiscard]] const WorkingMatrix& getMatrixU() const noexcept { assert(computeMatrixU); return matrixU; }
    private:
        template<class AnyMatrix>
        static size_t activeWindowLower(LValueMatrix<AnyMatrix>& __restrict mat, size_t upper);
        void splitOffTwoRows(size_t index);
        void francisQR(size_t lower, size_t sub_order);
        void specialHessenburg(size_t lower, size_t sub_order);
        ComplexType complexShift(size_t upper, size_t iter);
        void complexQR(size_t lower, size_t upper, ComplexType shift);
    };

    template<class MatrixType>
    template<class OtherMatrix>
    Schur<MatrixType>::Schur(const RValueMatrix<OtherMatrix>& source, bool computeMatrixU_)
            : matrixT(source.getRow(), source.getColumn())
            , matrixU()
            , computeMatrixU(computeMatrixU_) {
        assert(source.getRow() == source.getColumn());
        if (computeMatrixU)
            matrixU = WorkingMatrix::unitMatrix(source.getRow());

        typename MatrixType::RealMatrix buffer = abs(source);
        const RealType factor = buffer.max();
        if (factor < std::numeric_limits<ScalarType>::min()) {
            matrixT = RealType::Zero();
            return;
        }
        const RealType inv_factor = reciprocal(factor);
        const MatrixType normalized = source * inv_factor; //Referenced from eigen, guess to avoid overflow in householder, but will lost relative accuracy(from 10^-15 to 10^-14)
        const Hessenburg hess(normalized);
        matrixT = hess.getMatrixH();

        const size_t order = matrixT.getRow();
        size_t upper = order - 1;
        size_t iter = 0;
        size_t total_iter = 0;
        const size_t max_iter = matItePerCol * order;
        while (1 <= upper && upper < order) {
            const size_t lower = activeWindowLower(matrixT, upper);
            if (lower == upper) {
                upper -= 1;
                iter = 0;
            }
            else {
                if constexpr (ScalarType::isComplex) {
                    complexQR(lower, upper, complexShift(upper, iter));
                    ++iter;
                    ++total_iter;
                }
                else {
                    if (lower + 1 == upper) {
                        splitOffTwoRows(lower);
                        upper -= 2;
                        iter = 0;
                    }
                    else {
                        const size_t sub_order = upper - lower + 1;
                        francisQR(lower, sub_order);
                        ++iter;
                        ++total_iter;
                    }
                }
            }

            if (total_iter == max_iter)
                throw BadConvergenceException();
        }
        matrixT *= factor;

        if (computeMatrixU) {
            WorkingMatrix temp = WorkingMatrix(hess.getMatrixQ()) * matrixU;
            matrixU = std::move(temp);
        }
    }
    /**
     * We should process columns whose index is less or equal than \param upper
     * 
     * \returns We should process columns whose index is greater or equal to the return index
     */
    template<class MatrixType>
    template<class AnyMatrix>
    size_t Schur<MatrixType>::activeWindowLower(LValueMatrix<AnyMatrix>& __restrict mat, size_t upper) {
        assert(upper < mat.getRow());
        size_t lower = upper;
        size_t lower_1 = upper - 1;
        for (; lower_1 < lower; --lower, --lower_1) { //Make use of overflow
            RealType temp = abs(mat(lower, lower)) + abs(mat(lower_1, lower_1));
            temp = std::max(abs(temp * RealType(std::numeric_limits<ScalarType>::epsilon())), RealType(std::numeric_limits<ScalarType>::min()));
            if (abs(mat(lower, lower_1)) < temp) {
                mat(lower, lower_1) = ScalarType::Zero();
                break;
            }
        }
        return lower;
    }
    /**
     * Upper triangulize submatrix of \param mat, whose columns have index \param index and \param index + 1.
     */
    template<class MatrixType>
    void Schur<MatrixType>::splitOffTwoRows(size_t index) {
        const size_t index_1 = index + 1;
        const ScalarType p = ScalarType(0.5) * (matrixT(index, index) - matrixT(index_1, index_1));
        const ScalarType q = square(p) + matrixT(index, index_1) * matrixT(index_1, index);
        if (!q.isNegative()) {
            const ScalarType z = sqrt(q);
            Vector<ScalarType, 2> target;
            target[0] = p + (p.isPositive() ? z : -z); //Select the root that ensure numerical stable
            target[1] = matrixT(index_1, index);
            auto givensVector = givens(target, 0, 1);
            auto block1 = matrixT.rightCols(index);
            applyGivens(givensVector, block1, index, index_1);
            auto block2 = matrixT.topRows(index_1 + 1);
            givensVector[1].toOpposite();
            applyGivens(block2, givensVector, index, index_1);
            matrixT(index_1, index) = ScalarType::Zero();
            if (computeMatrixU)
                applyGivens(matrixU, givensVector, index, index_1);
        }
    }

    template<class MatrixType>
    void Schur<MatrixType>::francisQR(size_t lower, size_t sub_order) {
        auto subBlock = matrixT.block(lower, sub_order, lower, sub_order);
        const ScalarType s = subBlock(sub_order - 1, sub_order - 1) + subBlock(sub_order - 2, sub_order - 2);
        const ScalarType t = subBlock(sub_order - 1, sub_order - 1) * subBlock(sub_order - 2, sub_order - 2)
                           - subBlock(sub_order - 1, sub_order - 2) * subBlock(sub_order - 2, sub_order - 1);
        Vector<ScalarType, 3> col_1_M{};
        col_1_M[0] = square(subBlock(0, 0)) + subBlock(0, 1) * subBlock(1, 0) - s * subBlock(0, 0) + t;
        col_1_M[1] = subBlock(1, 0) * (subBlock(0, 0) + subBlock(1, 1) - s);

        if (sub_order != 2) {
            Vector<ScalarType, 3> householderVector{};
            col_1_M[2] = subBlock(1, 0) * subBlock(2, 1);
            householder(col_1_M, householderVector);
            {
                auto block = matrixT.rightCols(lower);
                auto rows = block.rows(lower, 3);
                applyHouseholder(householderVector, rows);
            }
            {
                auto block = matrixT.topRows(lower + sub_order);
                auto cols = block.cols(lower, 3);
                applyHouseholder(cols, householderVector);
            }
            if (computeMatrixU) {
                auto cols = matrixU.cols(lower, 3);
                applyHouseholder(cols, householderVector);
            }
            specialHessenburg(lower, sub_order);
        }
        else {
            Vector<ScalarType, 2> householderVector{};
            householder(col_1_M.head(2), householderVector);
            {
                auto block = matrixT.rightCols(lower);
                auto rows = block.rows(lower, 2);
                applyHouseholder(householderVector, rows);
            }
            {
                auto block = matrixT.topRows(lower + 2);
                auto cols = block.cols(lower, 2);
                applyHouseholder(cols, householderVector);
            }
            if (computeMatrixU) {
                auto cols = matrixU.cols(lower, 2);
                applyHouseholder(cols, householderVector);
            }
        }
    }
    /**
     * A special designed Hessenburg decomposition for francis QR algorithm
     */
    template<class MatrixType>
    void Schur<MatrixType>::specialHessenburg(size_t lower, size_t sub_order) {
        assert(sub_order > 2);
        Vector<ScalarType, 3> householderVector3D{};
        for (size_t i = 0; i < sub_order - 3; ++i) {
            auto block = matrixT.rows(lower + i + 1, 3);
            auto target_col = block.col(lower + i);
            const auto norm = householder(target_col, householderVector3D);
            target_col[0] = target_col[0].isNegative() ? norm : -norm;
            target_col.tail(1) = ScalarType::Zero();
            {
                auto rightCols = matrixT.rightCols(lower + i + 1);
                auto rows = rightCols.rows(lower + i + 1, 3);
                applyHouseholder(householderVector3D, rows);
            }
            {
                auto topRows = matrixT.topRows(lower + sub_order);
                auto cols = topRows.cols(lower + i + 1, 3);
                applyHouseholder(cols, householderVector3D);
            }
            if (computeMatrixU) {
                auto cols = matrixU.cols(lower + i + 1, 3);
                applyHouseholder(cols, householderVector3D);
            }
        }
        auto householderVector2D = householderVector3D.head(2);
        auto block = matrixT.rows(lower + sub_order - 2, 2);
        auto target_col = block.col(lower + sub_order - 3);
        const auto norm = householder(target_col, householderVector2D);
        target_col[0] = target_col[0].isNegative() ? norm : -norm;
        target_col[1] = ScalarType::Zero();
        {
            auto rightCols = matrixT.rightCols(lower + sub_order - 2);
            auto rows = rightCols.rows(lower + sub_order - 2, 2);
            applyHouseholder(householderVector2D, rows);
        }
        {
            auto topRows = matrixT.topRows(lower + sub_order);
            auto cols = topRows.cols(lower + sub_order - 2, 2);
            applyHouseholder(cols, householderVector2D);
        }
        if (computeMatrixU) {
            auto cols = matrixU.cols(lower + sub_order - 2, 2);
            applyHouseholder(cols, householderVector2D);
        }
    }

    template<class MatrixType>
    typename Schur<MatrixType>::ComplexType Schur<MatrixType>::complexShift(size_t upper, size_t iter) {
        using Matrix2D = DenseMatrix<ScalarType, DenseMatrixOption::Column | DenseMatrixOption::Element, 2, 2>;
        if ((iter == 10 || iter == 20) && upper > 1) {
            // exceptional shift, taken from http://www.netlib.org/eispack/comqr.f
            return abs(matrixT(upper, upper - 1).getReal()) + abs(matrixT(upper - 1, upper - 2).getReal());
        }
        //compute the shift as one of the eigenvalues of t, the 2x2
        //diagonal block on the bottom of the active submatrix
        const auto activeBlock = matrixT.block(upper - 1, 2, upper - 1, 2);
        RealType t_norm = abs(activeBlock(0, 0)) + abs(activeBlock(0, 1)) + abs(activeBlock(1, 0)) + abs(activeBlock(1, 1));
        const Matrix2D t = activeBlock * reciprocal(t_norm); //Normalization to avoid under/overflow

        const ComplexType b = t(0,1) * t(1,0);
        const ComplexType c = t(0,0) - t(1,1);
        const ComplexType disc = sqrt(square(c) + RealType(4) * b);
        const ComplexType det = t(0,0) * t(1,1) - b;
        const ComplexType trace = t(0,0) + t(1,1);
        ComplexType eival1 = (trace + disc) * RealType(0.5);
        ComplexType eival2 = (trace - disc) * RealType(0.5);
        const RealType eival1_norm = eival1.norm();
        const RealType eival2_norm = eival2.norm();

        if(eival1_norm > eival2_norm)
            eival2 = det / eival1;
        else if(!eival2_norm.isZero())
            eival1 = det / eival2;

        const bool firstEigenValueCloserToDiagonal = (eival1 - t(1,1)).norm() < (eival2 - t(1,1)).norm();
        return t_norm * (firstEigenValueCloserToDiagonal ? eival1 : eival2);
    }

    template<class MatrixType>
    void Schur<MatrixType>::complexQR(size_t lower, size_t upper, ComplexType shift) {
        using Vector2D = Vector<ScalarType, 2, 2>;
        {
            auto givensVec = givens(Vector2D{matrixT(lower, lower) - shift, matrixT(lower + 1, lower)}, 0, 1);
            auto rightCols = matrixT.rightCols(lower);
            applyGivens(givensVec, rightCols, lower, lower + 1);
            givensVec[1] = -givensVec[1];
            auto topRows = matrixT.topRows(std::min(lower + 2, upper) + 1);
            applyGivens(topRows, givensVec, lower, lower + 1);
            if (computeMatrixU)
                applyGivens(matrixU, givensVec, lower, lower + 1);
        }

        for(size_t i = lower + 1; i < upper; ++i) {
            auto givensVec = givens(Vector2D{matrixT(i, i - 1), matrixT(i + 1, i - 1)}, 0, 1);
            auto rightCols = matrixT.rightCols(i - 1);
            applyGivens(givensVec, rightCols, i, i + 1);
            matrixT(i + 1, i - 1) = ScalarType::Zero();
            givensVec[1] = -givensVec[1];
            auto topRows = matrixT.topRows(std::min(i + 2, upper) + 1);
            applyGivens(topRows, givensVec, i, i + 1);
            if (computeMatrixU)
                applyGivens(matrixU, givensVec, i, i + 1);
        }
    }
}
