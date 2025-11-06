/*
 * Copyright 2022-2025 Weibo He.
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

namespace Physica {
    template<Scalar T, int Option> class SparseMatrix;

    template<Scalar T, int Option, Vector V>
    class SpMV : public RValueVector<SpMV<T, Option, V>> {
        using Base = RValueVector<SpMV<T, Option, V>>;
    public:
        using typename Base::ScalarType;
        using MatrixType = SparseMatrix<T, Option>;
        static_assert(MatrixType::ColAtCompile == V::SizeAtCompile,
                      "Row and column do not match in matrix product");
    private:
        const MatrixType& mat;
        const V& vec;
    public:
        SpMV(const MatrixType& mat_, const V& vec_)
                : mat(mat_), vec(vec_) {
            assert(mat.getCol() == vec.getLength());
        }
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto& target) const;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const;
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow(); }
        [[nodiscard]] const auto& getLHS() const noexcept { return mat; }
        [[nodiscard]] const auto& getRHS() const noexcept { return vec; }
    };

    template<Scalar T, int Option, Vector V>
    template<ExecutePolicy P>
    void SpMV<T, Option, V>::assign(Vector auto& target) const {
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

    template<Scalar T, int Option, Vector V>
    auto SpMV<T, Option, V>::calc(size_t index) const -> ScalarType {
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
