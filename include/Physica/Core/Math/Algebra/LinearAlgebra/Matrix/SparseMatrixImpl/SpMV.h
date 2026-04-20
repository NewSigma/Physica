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

#include "../SparseMatrix.h"

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof_tx<SparseMatrix, M>)
    class GEMV<M, V> : public RValueVector<GEMV<M, V>> {
        using Base = RValueVector<GEMV<M, V>>;
    protected:
        using typename Base::T;
    private:
        LazyDestroy<M> mat;
        LazyDestroy<V> vec;
    public:
        GEMV(M&& mat_, V&& vec_);
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto& target) const;

        [[nodiscard]] T calc(size_t index) const;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow(); }
        [[nodiscard]] const auto& getLHS() const noexcept { return mat; }
        [[nodiscard]] const auto& getRHS() const noexcept { return vec; }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return Dynamic; }
    };

    template<Matrix M, Vector V> requires(instanceof_tx<SparseMatrix, M>)
    GEMV<M, V>::GEMV(M&& mat_, V&& vec_) : mat(std::forward<M>(mat_)), vec(std::forward<V>(vec_)) {
        assert(mat.getCol() == vec.getLength());
    }

    template<Matrix M, Vector V> requires(instanceof_tx<SparseMatrix, M>)
    template<ExecutePolicy P>
    void GEMV<M, V>::assign(Vector auto& target) const {
        const auto& elements = mat.getElements();
        const auto& minorIndexes = mat.getMinorIndexes();
        const auto& majorStarts = mat.getMajorStarts();

        if constexpr (MatrixMajor::isColMatrix<M>())
            target = 0;

        for (size_t major = 0; major < mat.getMaxMajor(); ++major) {
            const size_t from = majorStarts[major];
            const size_t to = majorStarts[major + 1];
            if constexpr (MatrixMajor::isColMatrix<M>()) {
                for (size_t i = from; i < to; ++i) {
                    const size_t row = minorIndexes[i];
                    target[row] += elements[i] * vec[major];
                }
            }
            else {
                T temp = 0;
                for (size_t i = from; i < to; ++i) {
                    const size_t col = minorIndexes[i];
                    temp += elements[i] * vec[col];
                }
                target[major] = temp;
            }
        }
    }

    template<Matrix M, Vector V> requires(instanceof_tx<SparseMatrix, M>)
    auto GEMV<M, V>::calc(size_t index) const -> T {
        const auto& elements = mat.getElements();
        const auto& minorIndexes = mat.getMinorIndexes();
        const auto& majorStarts = mat.getMajorStarts();

        T result = 0;
        if constexpr (MatrixMajor::isColMatrix<M>()) {
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
