/*
 * Copyright 2022 Weibo He.
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

#include "GEMV.h"

namespace Physica::Core {
    template<Scalar T, int Option> class SparseMatrix;

    template<Scalar T, int Option, Vector U>
    class MatrixVectorProduct<SparseMatrix<T, Option>, U>
            : public RValueVector<MatrixVectorProduct<SparseMatrix<T, Option>, U>> {
    public:
        using Base = RValueVector<MatrixVectorProduct<SparseMatrix<T, Option>, U>>;
        using typename Base::ScalarType;
        using MatrixType = SparseMatrix<T, Option>;
        static_assert(MatrixType::ColAtCompile == U::SizeAtCompile,
                      "Row and column do not match in matrix product");
    private:
        const MatrixType& mat;
        const U& vec;
    public:
        MatrixVectorProduct(const MatrixType& mat_, const U& vec_)
                : mat(mat_), vec(vec_) {
            assert(mat.getCol() == vec.getLength());
        }
        /* Operations */
        template<LVector V, class Executor = SequentialExecutor>
        void assignTo(V& target) const;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const;
        [[nodiscard]] size_t getLength() const { return mat.getRow(); }
        [[nodiscard]] const MatrixType& getLHS() const noexcept { return mat; }
        [[nodiscard]] const U& getRHS() const noexcept { return vec; }
    };

    template<Scalar T, int Option, Vector U>
    template<LVector V, class Executor>
    void MatrixVectorProduct<SparseMatrix<T, Option>, U>::assignTo(V& target) const {
        const auto& elements = mat.getElements();
        const auto& minorIndexes = mat.getMinorIndexes();
        const auto& majorStarts = mat.getMajorStarts();

        if constexpr (MatrixOption::isColMatrix<MatrixType>())
            target = 0;

        for (size_t major = 0; major < mat.getMaxMajor(); ++major) {
            const size_t from = majorStarts[major];
            const size_t to = majorStarts[major + 1];
            if constexpr (MatrixOption::isColMatrix<MatrixType>()) {
                for (size_t i = from; i < to; ++i) {
                    const size_t row = minorIndexes[i];
                    target[row] += elements[i] * vec[major];
                }
            }
            else {
                ScalarType temp = 0;
                for (size_t i = from; i < to; ++i) {
                    const size_t col = minorIndexes[i];
                    temp += elements[i] * vec[col];
                }
                target[major] = temp;
            }
        }
    }

    template<Scalar T, int Option, Vector U>
    typename MatrixVectorProduct<SparseMatrix<T, Option>, U>::ScalarType
    MatrixVectorProduct<SparseMatrix<T, Option>, U>::calc(size_t index) const {
        const auto& elements = mat.getElements();
        const auto& minorIndexes = mat.getMinorIndexes();
        const auto& majorStarts = mat.getMajorStarts();

        ScalarType result = 0;
        if constexpr (MatrixOption::isColMatrix<MatrixType>()) {
            for (size_t major = 0; major < mat.getMaxMajor(); ++major) {
                const size_t from = majorStarts[major];
                const size_t to = majorStarts[major + 1];
                for (size_t i = from; i < to; ++i) {
                    const size_t row = minorIndexes[i];
                    if (row == index)
                        result += elements[i] * vec[major];
                }
            }
        }
        else {
            const size_t from = majorStarts[index];
            const size_t to = majorStarts[index + 1];
            for (size_t i = from; i < to; ++i) {
                const size_t col = minorIndexes[i];
                result += elements[i] * vec[col];
            }
        }
        return result;
    }
}
