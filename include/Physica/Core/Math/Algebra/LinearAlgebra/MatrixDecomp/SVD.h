/*
 * Copyright 2022-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Givens.h"
#include "Physica/Core/Exception/BadConvergenceException.h"
#include "Bidiagonalization.h"
#include "Deflatable.h"

namespace Physica {
    /**
     * Decomposite matrix A like A = UDV^T
     *
     * References:
     * [1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013:488-492
     */
    template<Scalar T, size_t RowAtCompile = Dynamic, size_t ColAtCompile = Dynamic>
    class SVD : public Deflatable {
        static_assert(!T::isComplex(), "[Error]: SVD class do not support complex data");
        using Base = Deflatable;
        using Tr = T::RealType;
        constexpr static int Major = MatrixMajor::Col;
    public:
        constexpr static size_t NumSingularValue = RowAtCompile > ColAtCompile ? ColAtCompile : RowAtCompile;
        using SingularValueVector = DenseVector<Tr, NumSingularValue>;
        using WorkingMatrix = DenseMatrix<Tr, Major, RowAtCompile, ColAtCompile>;
        using LSingularMatrix = DenseMatrix<Tr, Major, RowAtCompile, RowAtCompile>;
        using RSingularMatrix = DenseMatrix<Tr, Major, ColAtCompile, ColAtCompile>;
    private:
        WorkingMatrix working;
        SingularValueVector singulars;
        LSingularMatrix lSingularMat;
        RSingularMatrix rSingularMat;
    public:
        SVD() = default;
        SVD(size_t row, size_t col);
        SVD(const Matrix auto& source);
        SVD(const SVD&) = default;
        SVD(SVD&&) noexcept = default;
        ~SVD() = default;
        /* Operators */
        SVD& operator=(SVD obj) noexcept;
        /* Operations */
        void compute(const Matrix auto& source);
        void sort();
        void resize(size_t row, size_t col);
        void resize(size_t order);
        void swap(SVD& __restrict svd) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getSignSingulars() const noexcept { return singulars; }
        [[nodiscard]] const auto& getMatrixU() const noexcept { return lSingularMat; }
        [[nodiscard]] const auto& getMatrixV() const noexcept { return rSingularMat; }
        [[nodiscard]] auto getSingulars() const noexcept { return abs(singulars); }
        [[nodiscard]] size_t getNumSingular() const noexcept { return singulars.getLength(); }
    private:
        void stepSVD(size_t lower, size_t sub_order);
        T computeShift(size_t lower, size_t sub_order);
        void leftGivens(Matrix auto& mat, Vector2D<T>& buffer, size_t i);
        void rightGivens(Matrix auto& mat, Vector2D<T>& buffer, size_t i);
    };

    template<Scalar T, size_t RowAtCompile, size_t ColAtCompile>
    SVD<T, RowAtCompile, ColAtCompile>::SVD(size_t row, size_t col) {
        resize(row, col);
    }

    template<Scalar T, size_t RowAtCompile, size_t ColAtCompile>
    SVD<T, RowAtCompile, ColAtCompile>::SVD(const Matrix auto& source)
            : SVD(source.getRow(), source.getCol()) {
        compute(source);
    }

    template<Scalar T, size_t RowAtCompile, size_t ColAtCompile>
    auto SVD<T, RowAtCompile, ColAtCompile>::operator=(SVD obj) noexcept -> SVD& {
        swap(obj);
        return *this;
    }

    template<Scalar T, size_t RowAtCompile, size_t ColAtCompile>
    void SVD<T, RowAtCompile, ColAtCompile>::compute(const Matrix auto& source) {
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
        const size_t max_iter = Deflatable::MaxIterationPerCol * order;
        const T factor = working.diag().normInf();
        while (1 <= upper && upper < order) {
            size_t lower = Base::activeWindowUpDiag(working, upper);
            for (size_t i = lower; i <= upper; ++i) {
                if (abs(working[i, i]) <= factor * T(std::numeric_limits<T>::epsilon())) {
                    working[i, i] = T(0);
                    working[i, i + (i + 1 < order)] = T(0);
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

            if (lower == upper)
                upper -= 1;
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
            singulars[i] = working[i, i];

        if (row < col)
            lSingularMat.swap(rSingularMat);
    }

    template<Scalar T, size_t RowAtCompile, size_t ColAtCompile>
    void SVD<T, RowAtCompile, ColAtCompile>::sort() {
        const size_t order = singulars.getLength();
        auto& larr = lSingularMat.asArray();
        auto& rarr = rSingularMat.asArray();
        for (size_t i = 0; i < order - 1; ++i) {
            const size_t index_min = i + abs(singulars.tail(i).reals()).argmin();
            if (index_min != i) {
                singulars[i].swap(singulars[index_min]);
                larr[i].swap(larr[index_min]);
                rarr[i].swap(rarr[index_min]);
            }
        }
    }

    template<Scalar T, size_t RowAtCompile, size_t ColAtCompile>
    void SVD<T, RowAtCompile, ColAtCompile>::resize(size_t row, size_t col) {
        working.resize(row, col);
        singulars.resize(std::min(row, col));
        lSingularMat.resize(row, row);
        rSingularMat.resize(col, col);
    }

    template<Scalar T, size_t RowAtCompile, size_t ColAtCompile>
    void SVD<T, RowAtCompile, ColAtCompile>::resize(size_t order) {
        resize(order, order);
    }

    template<Scalar T, size_t RowAtCompile, size_t ColAtCompile>
    void SVD<T, RowAtCompile, ColAtCompile>::swap(SVD& __restrict svd) noexcept {
        assert(this != &svd && "[Error]: Self swap is likely a bug");
        working.swap(svd.working);
        singulars.swap(svd.singulars);
        lSingularMat.swap(svd.lSingularMat);
        rSingularMat.swap(svd.rSingularMat);
    }

    template<Scalar T, size_t RowAtCompile, size_t ColAtCompile>
    void SVD<T, RowAtCompile, ColAtCompile>::stepSVD(size_t lower, size_t sub_order) {
        auto subBlock = working.block(lower, sub_order, lower, sub_order);
        const T shift = computeShift(lower, sub_order);

        Vector2D<T> buffer{square(subBlock[0, 0]) - shift, subBlock[0, 0] * subBlock[0, 1]};
        rightGivens(subBlock, buffer, 0);
        for (size_t i = 0; i < sub_order - 2; ++i) {
            leftGivens(subBlock, buffer, i);

            buffer[0] = subBlock[i, i + 1];
            buffer[1] = subBlock[i, i + 2];
            rightGivens(subBlock, buffer, i + 1);
            subBlock[i, i + 2] = 0;
        }
        leftGivens(subBlock, buffer, sub_order - 2);
    }
    /**
     * Use wilkinson shift
     */
    template<Scalar T, size_t RowAtCompile, size_t ColAtCompile>
    T SVD<T, RowAtCompile, ColAtCompile>::computeShift(size_t lower, size_t sub_order) {
        const size_t index = lower + sub_order - 1;
        const T d1 = working[index - 1, index - 1];
        const T d2 = working[index, index];
        T f1;
        if (lower == 0 && sub_order == 2)
            f1 = 0;
        else
            f1 = working[index - 2, index - 1];
        const T f2 = working[index - 1, index];
        const T a1 = square(d1) + square(f1);
        const T a2 = square(d2) + square(f2);
        const T b = d1 * f2;
        const Tr factor = (a1 - a2) * 0.5;
        const Tr factor2 = square(b);
        const Tr factor3 = sqrt(fma(factor, factor, factor2));
        const T shift = a2 - factor2 / (factor + (factor.isPositive() ? factor3 : -factor3));
        return shift;
    }

    template<Scalar T, size_t RowAtCompile, size_t ColAtCompile>
    void SVD<T, RowAtCompile, ColAtCompile>::leftGivens(Matrix auto& mat, Vector2D<T>& buffer, size_t i) {
        buffer[0] = mat[i, i];
        buffer[1] = mat[i + 1, i];
        buffer = givens(buffer, 0, 1);
        applyGivens(buffer, mat, i, i + 1);
        mat[i + 1, i] = 0;
        buffer[1] = -buffer[1];
        applyGivens(lSingularMat, buffer, i, i + 1);
    }

    template<Scalar T, size_t RowAtCompile, size_t ColAtCompile>
    void SVD<T, RowAtCompile, ColAtCompile>::rightGivens(Matrix auto& mat, Vector2D<T>& buffer, size_t i) {
        buffer = givens(buffer, 0, 1);
        buffer[1] = -buffer[1];
        applyGivens(mat, buffer, i, i + 1);
        applyGivens(rSingularMat, buffer, i, i + 1);
    }

    template<Scalar T, size_t RowAtCompile, size_t ColAtCompile>
    void swap(SVD<T, RowAtCompile, ColAtCompile>& svd1, SVD<T, RowAtCompile, ColAtCompile>& svd2) noexcept {
        svd1.swap(svd2);
    }
}
