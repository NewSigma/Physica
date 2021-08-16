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

#include "Hessenburg.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Givens.h"

namespace Physica::Core {
    template<class MatrixType> class RealSchur;

    namespace Internal {
        template<class MatrixType>
        class Traits<RealSchur<MatrixType>> : public Traits<MatrixType> {
        private:
            using Base = Traits<MatrixType>;
            using Base::MatrixOption;
        };
    }
    /**
     * References:
     * [1] Golub, GeneH. Matrix computations = 矩阵计算 / 4th edition[M]. 人民邮电出版社, 2014.
     * [2] Eigen https://eigen.tuxfamily.org/
     */
    template<class MatrixType>
    class RealSchur : public RValueMatrix<RealSchur<MatrixType>> {
        using Base = RValueMatrix<RealSchur<MatrixType>>;
        using typename Base::ScalarType;
        constexpr static size_t matItePerCol = 40; //Reference to Eigen
        const MatrixType& source;
    public:
        RealSchur(const LValueMatrix<MatrixType>& source_) : source(source_.getDerived()) { assert(source.getRow() == source.getColumn()); }
        /* Operations */
        template<class OtherMatrix>
        void assignTo(LValueMatrix<OtherMatrix>& target) const;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return source.getRow(); }
        [[nodiscard]] size_t getColumn() const noexcept { return source.getColumn(); }
    private:
        template<class AnyMatrix>
        static size_t activeWindowLower(LValueMatrix<AnyMatrix>& __restrict mat, size_t upper);
        template<class AnyMatrix>
        static void splitOffTwoRows(LValueMatrix<AnyMatrix>& __restrict mat, size_t index);
        template<class AnyMatrix>
        static void francisQR(LValueMatrix<AnyMatrix>& mat);
        template<class AnyMatrix>
        static void specialHessenburg(LValueMatrix<AnyMatrix>& mat);
    };

    template<class MatrixType>
    template<class OtherMatrix>
    void RealSchur<MatrixType>::assignTo(LValueMatrix<OtherMatrix>& target) const {
        Hessenburg hess(source);
        target = hess.getMatrixH();
        const size_t order = target.getRow();
        size_t upper = order - 1;
        size_t iter = 0;
        while (1 <= upper && upper < order) {
            const size_t lower = activeWindowLower(target, upper);
            if (lower == upper) {
                upper -= 1;
                iter = 0;
            }
            else if (lower + 1 == upper) {
                splitOffTwoRows(target, lower);
                upper -= 2;
                iter = 0;
            }
            else {
                const size_t sub_order = upper - lower + 1;
                auto subBlock = target.block(lower, sub_order, lower, sub_order);
                francisQR(subBlock);
                ++iter;
            }
        }
    }
    /**
     * We should process columns whose index is less or equal than \param upper
     * 
     * \returns We should process columns whose index is greater or equal to the return index
     */
    template<class MatrixType>
    template<class AnyMatrix>
    size_t RealSchur<MatrixType>::activeWindowLower(LValueMatrix<AnyMatrix>& __restrict mat, size_t upper) {
        assert(upper < mat.getRow());
        size_t lower = upper;
        size_t lower_1 = upper - 1;
        for (; lower_1 < lower; --lower, --lower_1) { //Make use of overflow
            ScalarType temp = abs(mat(lower, lower)) + abs(mat(lower_1, lower_1));
            temp = std::max(temp * std::numeric_limits<ScalarType>::epsilon(), ScalarType(std::numeric_limits<ScalarType>::min()));
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
    template<class AnyMatrix>
    void RealSchur<MatrixType>::splitOffTwoRows(LValueMatrix<AnyMatrix>& __restrict mat, size_t index) {
        const size_t index_1 = index + 1;
        const ScalarType p = ScalarType(0.5) * (mat(index, index) - mat(index_1, index_1));
        const ScalarType q = square(p) + mat(index, index_1) * mat(index_1, index);
        if (!q.isNegative()) {
            const ScalarType z = sqrt(q);
            Vector<ScalarType, 2> target;
            target[0] = p + (p.isPositive() ? z : -z); //Select the root that ensure numerical stable
            target[1] = mat(index_1, index);
            auto givensVector = givens(target, 0, 1);
            auto block1 = mat.rightCols(index);
            applyGivens(givensVector, block1, index, index_1);
            auto block2 = mat.topRows(index_1 + 1);
            givensVector[1].toOpposite();
            applyGivens(block2, givensVector, index, index_1);
            mat(index_1, index) = ScalarType::Zero();
        }
    }

    template<class MatrixType>
    template<class AnyMatrix>
    void RealSchur<MatrixType>::francisQR(LValueMatrix<AnyMatrix>& mat) {
        constexpr size_t householderSize = AnyMatrix::RowAtCompile == 2 ? 2 : 3;
        const size_t order = mat.getRow();
        const ScalarType s = mat(order - 1, order - 1) + mat(order - 2, order - 2);
        const ScalarType t = mat(order - 1, order - 1) * mat(order - 2, order - 2)
                           - mat(order - 1, order - 2) * mat(order - 2, order - 1);
        Vector<ScalarType, 3> col_1_M{};
        col_1_M[0] = square(mat(0, 0)) + mat(0, 1) * mat(1, 0) - s * mat(0, 0) + t;
        col_1_M[1] = mat(1, 0) * (mat(0, 0) + mat(1, 1) - s);

        Vector<ScalarType, householderSize> householderVector{};
        if (order != 2) {
            col_1_M[2] = mat(1, 0) * mat(2, 1);
            householder(col_1_M, householderVector);
            auto rows = mat.rows(0, 3);
            applyHouseholder(householderVector, rows);
            auto cols = mat.cols(0, 3);
            applyHouseholder(cols, householderVector);

            specialHessenburg(mat);
        }
        else {
            auto householderVector2D = householderVector.head(2);
            householder(col_1_M.head(2), householderVector2D);
            applyHouseholder(householderVector2D, mat);
            applyHouseholder(mat, householderVector2D);
        }
    }
    /**
     * A special designed Hessenburg decomposition for francis QR algorithm
     */
    template<class MatrixType>
    template<class AnyMatrix>
    void RealSchur<MatrixType>::specialHessenburg(LValueMatrix<AnyMatrix>& mat) {
        assert(mat.getRow() > 2);
        Vector<ScalarType, 3> householderVector3D{};
        const size_t order = mat.getRow();
        for (size_t i = 0; i < order - 3; ++i) {
            auto block = mat.rows(i + 1, 3);
            auto target_col = block.col(i);
            const auto norm = householder(target_col, householderVector3D);
            target_col[0] = target_col[0].isNegative() ? norm : -norm;
            target_col.tail(1) = ScalarType::Zero();
            auto rightCols = block.rightCols(i + 1);
            applyHouseholder(householderVector3D, rightCols);
            auto cols = mat.cols(i + 1, 3);
            applyHouseholder(cols, householderVector3D);
        }
        auto householderVector2D = householderVector3D.head(2);
        auto block = mat.bottomRows(order - 2);
        auto target_col = block.col(order - 3);
        const auto norm = householder(target_col, householderVector2D);
        target_col[0] = target_col[0].isNegative() ? norm : -norm;
        target_col[1] = ScalarType::Zero();
        auto rightCols = block.rightCols(order - 2);
        applyHouseholder(householderVector2D, rightCols);
        auto cols = mat.rightCols(order - 2);
        applyHouseholder(cols, householderVector2D);
    }
}
