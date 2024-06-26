/*
 * Copyright 2022 WeiBo He.
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

#include "MatrixProduct.h"

namespace Physica::Core {
    template<class T, int Option, class VectorType>
    class MatrixVectorProduct<SparseMatrix<T, Option>, VectorType>
            : public RValueVector<MatrixVectorProduct<SparseMatrix<T, Option>, VectorType>> {
    public:
        using Base = RValueVector<MatrixVectorProduct<SparseMatrix<T, Option>, VectorType>>;
        using typename Base::ScalarType;
        using MatrixType = SparseMatrix<T, Option>;
        static_assert(MatrixType::ColumnAtCompile == VectorType::SizeAtCompile,
                      "Row and column do not match in matrix product");
    private:
        const MatrixType& mat;
        const VectorType& vec;
    public:
        MatrixVectorProduct(const MatrixType& mat_, const RValueVector<VectorType>& vec_)
                : mat(mat_), vec(vec_.getDerived()) {
            assert(mat.getColumn() == vec.getLength());
        }
        /* Operations */
        template<class OtherDerived, class Executor = SequentialExecutor>
        void assignTo(LValueVector<OtherDerived>& target) const;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const;
        [[nodiscard]] size_t getLength() const { return mat.getRow(); }
        [[nodiscard]] const MatrixType& getLHS() const noexcept { return mat; }
        [[nodiscard]] const VectorType& getRHS() const noexcept { return vec; }
    };

    template<class T, int Option, class VectorType>
    template<class OtherDerived, class Executor>
    void MatrixVectorProduct<SparseMatrix<T, Option>, VectorType>::assignTo(LValueVector<OtherDerived>& target) const {
        const auto& elements = mat.getElements();
        const auto& minorIndexes = mat.getMinorIndexes();
        const auto& majorStarts = mat.getMajorStarts();

        if constexpr (MatrixOption::isColumnMatrix<MatrixType>())
            target = 0;

        for (size_t major = 0; major < mat.getMaxMajor(); ++major) {
            const size_t from = majorStarts[major];
            const size_t to = majorStarts[major + 1];
            if constexpr (MatrixOption::isColumnMatrix<MatrixType>()) {
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

    template<class T, int Option, class VectorType>
    typename MatrixVectorProduct<SparseMatrix<T, Option>, VectorType>::ScalarType
    MatrixVectorProduct<SparseMatrix<T, Option>, VectorType>::calc(size_t index) const {
        const auto& elements = mat.getElements();
        const auto& minorIndexes = mat.getMinorIndexes();
        const auto& majorStarts = mat.getMajorStarts();

        ScalarType result = 0;
        if constexpr (MatrixOption::isColumnMatrix<MatrixType>()) {
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
