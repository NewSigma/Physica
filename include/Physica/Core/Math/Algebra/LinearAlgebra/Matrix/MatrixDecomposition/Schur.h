/*
 * Copyright 2021-2024 Weibo He.
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
#include "Decouplable.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Givens.h"
#include "Physica/Core/Exception/BadConvergenceException.h"

namespace Physica::Core {
    /**
     * A = UTU^H
     * 
     * References:
     * [1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013
     * [2] Eigen; https://gitlab.com/libeigen/eigen
     */
    template<class ScalarType, size_t Order = Dynamic>
    class Schur : public Decouplable {
        static_assert(is_scalar<ScalarType>::value, "[Error]: Invalid ScalarType");
        constexpr static const char* BadConvergenceMessage = "Exceed max iteration of Schur";
        using HessenburgType = Hessenburg<ScalarType, Order>;
        using Vector3D = Vector<ScalarType, 3>;
    public:
        using RealType = typename ScalarType::RealType;
        using ComplexType = ComplexScalar<RealType>;
        using WorkingMatrix = typename HessenburgType::WorkingMatrix;
    private:
        WorkingMatrix matrixT;
        WorkingMatrix matrixU;
        bool computeMatrixU;
        ScalarType exshift;
    public:
        template<class MatrixType>
        Schur(const RValueMatrix<MatrixType>& source_, bool computeMatrixU_ = false);
        /* Getters */
        [[nodiscard]] WorkingMatrix& getMatrixT() noexcept { return matrixT; }
        [[nodiscard]] WorkingMatrix& getMatrixU() noexcept { assert(computeMatrixU); return matrixU; }
        [[nodiscard]] const WorkingMatrix& getMatrixT() const noexcept { return matrixT; }
        [[nodiscard]] const WorkingMatrix& getMatrixU() const noexcept { assert(computeMatrixU); return matrixU; }
    private:
        void splitOffTwoRows(size_t index);
        Vector3D realShift(size_t upper, size_t iter);
        void francisQR(size_t lower, size_t sub_order, Vector3D shift);
        void specialHessenburg(size_t lower, size_t sub_order);

        ComplexType complexShift(size_t upper, size_t iter);
        void complexQR(size_t lower, size_t upper, ComplexType shift);
    };

    template<class ScalarType, size_t Order>
    template<class MatrixType>
    Schur<ScalarType, Order>::Schur(const RValueMatrix<MatrixType>& source_, bool computeMatrixU_)
            : matrixT(source_.getRow(), source_.getColumn())
            , matrixU()
            , computeMatrixU(computeMatrixU_) {
        static_assert(std::is_same<ScalarType, typename MatrixType::ScalarType>::value, "[Error]: Inconsistent ScalarType");
        const auto& source = source_.getDerived();
        assert(source.getRow() == source.getColumn());
        if (computeMatrixU)
            matrixU = WorkingMatrix::unitMatrix(source.getRow());

        const RealType factor = abs_elem(source).max();
        if (factor < std::numeric_limits<ScalarType>::min()) {
            matrixT = RealType(0);
            return;
        }
        const RealType inv_factor = reciprocal(factor);
        const auto normalized = source * inv_factor; //Referenced from eigen, to avoid under/overflow in householder, but will lost relative accuracy(from 10^-15 to 10^-14)

        const size_t order = source.getRow();
        size_t iter = 0;
        exshift = 0;
        if (order != 2) {
            const Hessenburg<ScalarType, Order> hess(normalized);
            matrixT = hess.getMatrixH();

            size_t upper = order - 1;
            size_t total_iter = 0;
            const size_t max_iter = Decouplable::maxItePerCol * order;
            while (1 <= upper && upper < order) {
                const size_t lower = activeWindowDownDiag(matrixT, upper);
                if (lower == upper) {
                    matrixT(upper, upper) += exshift;
                    upper -= 1;
                    iter = 0;
                }
                else {
                    if (total_iter == max_iter) [[unlikely]]
                        throw BadConvergenceException(BadConvergenceMessage);

                    if constexpr (ScalarType::isComplex) {
                        const auto shift = complexShift(upper, iter);
                        complexQR(lower, upper, shift);
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
                            const auto shift = realShift(upper, iter);
                            francisQR(lower, sub_order, shift);
                            ++iter;
                            ++total_iter;
                        }
                    }
                }
            }
            matrixT(0, 0) += exshift;

            if (computeMatrixU) {
                WorkingMatrix temp = WorkingMatrix(hess.getMatrixQ()) * matrixU;
                matrixU = std::move(temp);
            }
        }
        else {
            matrixT = normalized;
            if constexpr (ScalarType::isComplex) {
                while (activeWindowDownDiag(matrixT, 1) != 1) {
                    if (iter == Decouplable::maxItePerCol) [[unlikely]]
                        throw BadConvergenceException(BadConvergenceMessage);
                    const auto shift = complexShift(1, iter);
                    complexQR(0, 1, shift);
                    iter += 1;
                }
            }
            else {
                if (activeWindowDownDiag(matrixT, 1) != 1)
                    splitOffTwoRows(0);
            }
        }
        matrixT *= factor;
    }
    /**
     * Upper triangulize submatrix of \param mat, whose columns have index \param index and \param index + 1.
     */
    template<class ScalarType, size_t Order>
    void Schur<ScalarType, Order>::splitOffTwoRows(size_t index) {
        const size_t index_1 = index + 1;
        matrixT(index, index) += exshift;
        matrixT(index_1, index_1) += exshift;

        const ScalarType p = ScalarType(0.5) * (matrixT(index, index) - matrixT(index_1, index_1));
        const ScalarType q = square(p) + matrixT(index, index_1) * matrixT(index_1, index);
        const bool haveTwoRealEigenvalues = !q.isNegative();
        if (haveTwoRealEigenvalues) {
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
            matrixT(index_1, index) = ScalarType(0);
            if (computeMatrixU)
                applyGivens(matrixU, givensVector, index, index_1);
        }
    }
    /**
     * Reference:
     * [1] Eigen; https://gitlab.com/libeigen/eigen
     */
    template<class ScalarType, size_t Order>
    typename Schur<ScalarType, Order>::Vector3D Schur<ScalarType, Order>::realShift(size_t upper, size_t iter) {
        const size_t upper_1 = upper - 1;
        ScalarType s = matrixT(upper, upper) + matrixT(upper_1, upper_1);
        ScalarType t1 = matrixT(upper, upper) * matrixT(upper_1, upper_1);
        ScalarType t2 = matrixT(upper, upper_1) * matrixT(upper_1, upper);
        if (iter > 0 && iter % 16 == 0) {
            const bool useWilkinsonShift = iter % 32 != 0;
            if (useWilkinsonShift) {
                const ScalarType shift = matrixT(upper, upper);
                exshift += shift;
                for (size_t i = 0; i <= upper; ++i)
                    matrixT(i, i) -= shift;

                const ScalarType s1 = abs(matrixT(upper, upper_1)) + abs(matrixT(upper_1, upper - 2));
                const ScalarType s2 = square(s1);
                s = ScalarType(0.75 + 0.75) * s1;
                t1 = ScalarType(0.75 * 0.75) * s2;
                t2 = ScalarType(-0.4375) * s2;
            }
            else { // MATLAB new ad hoc shift
                const ScalarType s1 = (matrixT(upper_1, upper_1) - matrixT(upper, upper)) * ScalarType(0.5);
                ScalarType shift = square(s1) + t2;
                if (s.isPositive()) {
                    shift = sqrt(shift);
                    if (s1.isNegative())
                        shift = -shift;
                    shift += s1;
                    shift = matrixT(upper, upper) - t2 / shift;
                    exshift += shift;
                    for (size_t i = 0; i <= upper; ++i)
                        matrixT(i, i) -= shift;
                    s = ScalarType(0.964 * 2);
                    t1 = ScalarType(0.964 * 0.964);
                    t2 = ScalarType(0.964);
                }
            }
        }
        return {s, t1, t2};
    }

    template<class ScalarType, size_t Order>
    void Schur<ScalarType, Order>::francisQR(size_t lower, size_t sub_order, Vector3D shift) {
        auto subBlock = matrixT.block(lower, sub_order, lower, sub_order);
        Vector3D col_1_M{};
        col_1_M[0] = (subBlock(0, 0) - shift[0]) * subBlock(0, 0) + shift[1] + (subBlock(0, 1) * subBlock(1, 0) - shift[2]);
        col_1_M[1] = subBlock(1, 0) * (subBlock(0, 0) + subBlock(1, 1) - shift[0]);

        if (sub_order != 2) {
            Vector3D householderVector{};
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
    template<class ScalarType, size_t Order>
    void Schur<ScalarType, Order>::specialHessenburg(size_t lower, size_t sub_order) {
        assert(sub_order > 2);
        Vector3D householderVector3D{};
        for (size_t i = 0; i < sub_order - 3; ++i) {
            auto block = matrixT.rows(lower + i + 1, 3);
            auto target_col = block.col(lower + i);
            const auto norm = householder(target_col, householderVector3D);
            target_col[0] = target_col[0].isNegative() ? norm : -norm;
            target_col.tail(1) = ScalarType(0);
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
        target_col[1] = ScalarType(0);
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

    template<class ScalarType, size_t Order>
    typename Schur<ScalarType, Order>::ComplexType Schur<ScalarType, Order>::complexShift(size_t upper, size_t iter) {
        using Matrix2D = DenseMatrix<ScalarType, MatrixOption::Column | MatrixOption::Element, 2, 2>;
        if ((iter == 10 || iter == 20) && upper > 1) {
            // exceptional shift, taken from http://www.netlib.org/eispack/comqr.f
            return abs(matrixT(upper, upper - 1).getReal()) + abs(matrixT(upper - 1, upper - 2).getReal());
        }
        //compute the shift as one of the eigenvalues of t, the 2x2
        //diagonal block on the bottom of the active submatrix
        const auto activeBlock = matrixT.block(upper - 1, 2, upper - 1, 2);
        const RealType t_norm = abs(activeBlock(0, 0)) + abs(activeBlock(0, 1)) + abs(activeBlock(1, 0)) + abs(activeBlock(1, 1));
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

    template<class ScalarType, size_t Order>
    void Schur<ScalarType, Order>::complexQR(size_t lower, size_t upper, ComplexType shift) {
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
            matrixT(i + 1, i - 1) = ScalarType(0);
            givensVec[1] = -givensVec[1];
            auto topRows = matrixT.topRows(std::min(i + 2, upper) + 1);
            applyGivens(topRows, givensVec, i, i + 1);
            if (computeMatrixU)
                applyGivens(matrixU, givensVec, i, i + 1);
        }
    }
}
