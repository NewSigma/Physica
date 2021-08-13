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
        static void updateActiveWindow(LValueMatrix<AnyMatrix>& __restrict mat, size_t& __restrict lower, size_t& __restrict upper);
        template<class AnyMatrix>
        static void francisQR(LValueMatrix<AnyMatrix>& mat);
        template<class AnyMatrix>
        static void specialHessenburg(LValueMatrix<AnyMatrix>& mat);
    };

    template<class MatrixType>
    template<class OtherMatrix>
    void RealSchur<MatrixType>::assignTo(LValueMatrix<OtherMatrix>& target) const {
        target = Hessenburg(source);
        const size_t order = target.getRow();
        //Upper and lower index for active window
        size_t lower_index = 0;
        size_t upper_index = order - 1;
        do {
            updateActiveWindow(target, lower_index, upper_index);
            francisQR(target);
        } while(lower_index != upper_index);
    }
    /**
     * We are processing columns whose index is greater or equal to \param lower and less than \param upper
     */
    template<class MatrixType>
    template<class AnyMatrix>
    void RealSchur<MatrixType>::updateActiveWindow(LValueMatrix<AnyMatrix>& __restrict mat, size_t& __restrict lower, size_t& __restrict upper) {
        for (size_t lower_1 = lower + 1; lower < upper; ++lower_1, ++lower) {
            ScalarType temp = abs(mat(lower, lower)) + abs(mat(lower_1, lower_1));
            temp = std::max(temp * std::numeric_limits<ScalarType>::epsilon(), ScalarType(std::numeric_limits<ScalarType>::min()));
            if (abs(mat(lower_1, lower)) >= temp)
                break;
            mat(lower_1, lower) = ScalarType::Zero();
        }

        for (size_t upper_1 = upper - 1; upper > lower; --upper_1, --upper) {
            ScalarType temp = abs(mat(upper_1, upper_1)) + abs(mat(upper, upper));
            temp = std::max(temp * std::numeric_limits<ScalarType>::epsilon(), ScalarType(std::numeric_limits<ScalarType>::min()));
            if (abs(mat(upper, upper_1)) >= temp)
                break;
            mat(upper, upper_1) = ScalarType::Zero();
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
        col_1_M[2] = mat(1, 0) * mat(2, 1);

        Vector<ScalarType, householderSize> householderVector{};
        if (order != 2)
            householder(col_1_M, householderVector);
        else
            householder(col_1_M.head(2), householderVector);
        auto rows = mat.rows(0, 3);
        applyHouseholder(householderVector, rows);
        auto cols = mat.cols(0, 3);
        applyHouseholder(cols, householderVector);

        specialHessenburg(mat);
    }
    /**
     * A special designed Hessenburg decomposition for francis QR algorithm
     */
    template<class MatrixType>
    template<class AnyMatrix>
    void RealSchur<MatrixType>::specialHessenburg(LValueMatrix<AnyMatrix>& mat) {
        Vector<ScalarType, 3> householderVector3D{};
        const size_t order = mat.getRow();
        if (order > 2) {
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
