/*
 * Copyright 2022-2024 Weibo He.
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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector/Givens.h>
#include <Physica/Core/Exception/BadConvergenceException.h>
#include "Bidiagonalization.h"
#include "Decouplable.h"

namespace Physica::Core {
    /**
     * Decomposite matrix A like A = UDV^T
     * 
     * References:
     * [1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013:488-492
     */
    template<class ScalarType, size_t RowAtCompile = Dynamic, size_t ColumnAtCompile = Dynamic>
    class SVD : public Decouplable {
        using Base = Decouplable;
        using RealType = typename ScalarType::RealType;
        static_assert(!ScalarType::isComplex, "[Error]: SVD class do not support complex data");
        constexpr static size_t NumSingularValue = RowAtCompile > ColumnAtCompile
                                                                ? ColumnAtCompile
                                                                : RowAtCompile;
        using WorkingMatrix = DenseMatrix<RealType,
                                          MatrixOption::Col | MatrixOption::Vector,
                                          RowAtCompile,
                                          ColumnAtCompile>;
    public:
        using SingularValueVector = Vector<RealType, NumSingularValue>;
        using LSingularMatrix = DenseMatrix<RealType,
                                            MatrixOption::Col | MatrixOption::Vector,
                                            RowAtCompile,
                                            RowAtCompile>;
        using RSingularMatrix = DenseMatrix<RealType,
                                            MatrixOption::Col | MatrixOption::Vector,
                                            ColumnAtCompile,
                                            ColumnAtCompile>;
    private:
        WorkingMatrix working;
        SingularValueVector singulars;
        LSingularMatrix lSingularMat;
        RSingularMatrix rSingularMat;
    public:
        SVD() = default;
        SVD(size_t row, size_t col);
        template<class OtherMatrix>
        SVD(const RValueMatrix<OtherMatrix>& source);
        SVD(const SVD&) = default;
        SVD(SVD&&) noexcept = default;
        ~SVD() = default;
        /* Operators */
        SVD& operator=(SVD obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class OtherMatrix>
        void compute(const RValueMatrix<OtherMatrix>& source);
        void sort();
        void resize(size_t row, size_t col);
        void swap(SVD& __restrict svd) noexcept;
        /* Getters */
        [[nodiscard]] const SingularValueVector& getSingulars() const noexcept { return singulars; }
        [[nodiscard]] const LSingularMatrix& getMatrixU() const noexcept { return lSingularMat; }
        [[nodiscard]] const RSingularMatrix& getMatrixV() const noexcept { return rSingularMat; }
        [[nodiscard]] size_t getNumSingular() const noexcept { return singulars.getLength(); }
    private:
        void stepSVD(size_t lower, size_t sub_order);
        ScalarType computeShift(size_t lower, size_t sub_order);
        template<class MatrixType>
        void leftGivens(LValueMatrix<MatrixType>& mat, Vector<ScalarType, 2>& buffer, size_t i);
        template<class MatrixType>
        void rightGivens(LValueMatrix<MatrixType>& mat, Vector<ScalarType, 2>& buffer, size_t i);
    };

    template<class ScalarType, size_t RowAtCompile, size_t ColumnAtCompile>
    SVD<ScalarType, RowAtCompile, ColumnAtCompile>::SVD(size_t row, size_t col) {
        resize(row, col);
    }

    template<class ScalarType, size_t RowAtCompile, size_t ColumnAtCompile>
    template<class OtherMatrix>
    SVD<ScalarType, RowAtCompile, ColumnAtCompile>::SVD(const RValueMatrix<OtherMatrix>& source)
            : SVD(source.getRow(), source.getCol()) {
        compute(source);
    }

    template<class ScalarType, size_t RowAtCompile, size_t ColumnAtCompile>
    template<class OtherMatrix>
    void SVD<ScalarType, RowAtCompile, ColumnAtCompile>::compute(const RValueMatrix<OtherMatrix>& source) {
        const size_t row = source.getRow();
        const size_t col = source.getCol();
        assert(row == lSingularMat.getRow());
        assert(col == rSingularMat.getRow());
        {
            Bidiagonalization<WorkingMatrix> bidiag(std::max(row, col), std::min(row, col));
            if (row < col) {
                bidiag.compute(source.transpose());
                rSingularMat = bidiag.getMatrixU();
                lSingularMat = bidiag.getMatrixV();
            }
            else {
                bidiag.compute(source);
                lSingularMat = bidiag.getMatrixU();
                rSingularMat = bidiag.getMatrixV();
            }
            working = bidiag.getMatrixB();
        }
        const size_t order = working.getCol();
        size_t upper = order - 1;
        size_t total_iter = 0;
        const size_t max_iter = Decouplable::maxItePerCol * order;
        while (1 <= upper && upper < order) {
            size_t lower = Base::activeWindowUpDiag(working, upper);
            for (size_t i = lower; i <= upper; ++i) {
                if (abs(working(i, i)) <= ScalarType(std::numeric_limits<ScalarType>::epsilon() * 10)) { //TODO: absolute criteria is unrelirable, use relative criteria instead 
                    working(i, i) = ScalarType(0);
                    working(i, i + (i + 1 < order)) = ScalarType(0);
                    if (i < upper - 1) {
                        lower = i + 1;
                        break;
                    }
                    if (i == upper - 1) {
                        upper -= 1;
                        goto pass;
                    }
                }
            }

            if (lower == upper) {
                upper -= 1;
            }
            else {
                const size_t sub_order = upper - lower + 1;
                stepSVD(lower, sub_order);
                ++total_iter;
            }

            if (total_iter == max_iter)
                throw BadConvergenceException("Exceed max iteration of SVD");
            pass:;
        }

        for (size_t i = 0; i < order; ++i)
            singulars[i] = working(i, i);
        if (row < col)
            lSingularMat.swap(rSingularMat);
    }

    template<class ScalarType, size_t RowAtCompile, size_t ColumnAtCompile>
    void SVD<ScalarType, RowAtCompile, ColumnAtCompile>::sort() {
        const size_t order = singulars.getLength();
        for (size_t i = 0; i < order - 1; ++i) {
            size_t index_min = i;
            for (size_t j = i + 1; j < order; ++j) {
                if (abs(singulars[j].real()) < abs(singulars[index_min].real()))
                    index_min = j;
            }
            singulars[i].swap(singulars[index_min]);
            lSingularMat[i].swap(lSingularMat[index_min]);
            rSingularMat[i].swap(rSingularMat[index_min]);
        }
    }

    template<class ScalarType, size_t RowAtCompile, size_t ColumnAtCompile>
    void SVD<ScalarType, RowAtCompile, ColumnAtCompile>::resize(size_t row, size_t col) {
        working.resize(row, col);
        singulars.resize(std::min(row, col));
        lSingularMat.resize(row, row);
        rSingularMat.resize(col, col);
    }

    template<class ScalarType, size_t RowAtCompile, size_t ColumnAtCompile>
    void SVD<ScalarType, RowAtCompile, ColumnAtCompile>::swap(SVD<ScalarType, RowAtCompile, ColumnAtCompile>& __restrict svd) noexcept {
        assert(this != &svd && "[Error]: Self swap is likely a bug");
        working.swap(svd.working);
        singulars.swap(svd.singulars);
        lSingularMat.swap(svd.lSingularMat);
        rSingularMat.swap(svd.rSingularMat);
    }

    template<class ScalarType, size_t RowAtCompile, size_t ColumnAtCompile>
    void SVD<ScalarType, RowAtCompile, ColumnAtCompile>::stepSVD(size_t lower, size_t sub_order) {
        auto subBlock = working.block(lower, sub_order, lower, sub_order);
        const ScalarType shift = computeShift(lower, sub_order);

        Vector<ScalarType, 2> buffer{square(subBlock(0, 0)) - shift, subBlock(0, 0) * subBlock(0, 1)};
        rightGivens(subBlock, buffer, 0);
        for (size_t i = 0; i < sub_order - 2; ++i) {
            leftGivens(subBlock, buffer, i);

            buffer[0] = subBlock(i, i + 1);
            buffer[1] = subBlock(i, i + 2);
            rightGivens(subBlock, buffer, i + 1);
            subBlock(i, i + 2) = 0;
        }
        leftGivens(subBlock, buffer, sub_order - 2);
    }
    /**
     * Use wilkinson shift
     */
    template<class ScalarType, size_t RowAtCompile, size_t ColumnAtCompile>
    ScalarType SVD<ScalarType, RowAtCompile, ColumnAtCompile>::computeShift(size_t lower, size_t sub_order) {
        const size_t index = lower + sub_order - 1;
        const ScalarType d1 = working(index - 1, index - 1);
        const ScalarType d2 = working(index, index);
        ScalarType f1;
        if (lower == 0 && sub_order == 2)
            f1 = 0;
        else
            f1 = working(index - 2, index - 1);
        const ScalarType f2 = working(index - 1, index);
        const ScalarType a1 = square(d1) + square(f1);
        const ScalarType a2 = square(d2) + square(f2);
        const ScalarType b = d1 * f2;
        const RealType factor = (a1 - a2) * 0.5;
        const RealType factor2 = square(b);
        const RealType factor3 = sqrt(square(factor) + factor2);
        const ScalarType shift = a2 - factor2 / (factor + (factor.isPositive() ? factor3 : -factor3));
        return shift;
    }

    template<class ScalarType, size_t RowAtCompile, size_t ColumnAtCompile>
    template<class MatrixType>
    void SVD<ScalarType, RowAtCompile, ColumnAtCompile>::leftGivens(
            LValueMatrix<MatrixType>& mat,
            Vector<ScalarType, 2>& buffer,
            size_t i) {
        buffer[0] = mat(i, i);
        buffer[1] = mat(i + 1, i);
        buffer = givens(buffer, 0, 1);
        applyGivens(buffer, mat, i, i + 1);
        mat(i + 1, i) = 0;
        buffer[1] = -buffer[1];
        applyGivens(lSingularMat, buffer, i, i + 1);
    }

    template<class ScalarType, size_t RowAtCompile, size_t ColumnAtCompile>
    template<class MatrixType>
    void SVD<ScalarType, RowAtCompile, ColumnAtCompile>::rightGivens(
            LValueMatrix<MatrixType>& mat,
            Vector<ScalarType, 2>& buffer,
            size_t i) {
        buffer = givens(buffer, 0, 1);
        buffer[1] = -buffer[1];
        applyGivens(mat, buffer, i, i + 1);
        applyGivens(rSingularMat, buffer, i, i + 1);
    }
}

namespace std {
    template<class ScalarType, size_t RowAtCompile, size_t ColumnAtCompile>
    inline void swap(Physica::Core::SVD<ScalarType, RowAtCompile, ColumnAtCompile>& __restrict svd1,
                     Physica::Core::SVD<ScalarType, RowAtCompile, ColumnAtCompile>& __restrict svd2) noexcept {
        svd1.swap(svd2);
    }
}
