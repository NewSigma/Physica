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

#include "RValueMatrix.h"
#include "MatrixOption.h"

namespace Physica::Core {
    template<class ScalarType, int option = MatrixOption::Row | MatrixOption::Element>
    class SparseMatrix;

    namespace Internal {
        template<class T, int option>
        class Traits<SparseMatrix<T, option>> {
        public:
            using ScalarType = T;
            constexpr static int Option = option;
            constexpr static size_t RowAtCompile = Utils::Dynamic;
            constexpr static size_t ColumnAtCompile = Utils::Dynamic;
            constexpr static size_t MaxRowAtCompile = Utils::Dynamic;
            constexpr static size_t MaxColumnAtCompile = Utils::Dynamic;
            constexpr static size_t SizeAtCompile = Utils::Dynamic;
            constexpr static size_t MaxSizeAtCompile = Utils::Dynamic;
        };
    }

    template<class ScalarType, int option>
    class SparseMatrix : public RValueMatrix<SparseMatrix<ScalarType, option>> {
        using This = SparseMatrix<ScalarType, option>;
        using Base = RValueMatrix<This>;
        static_assert(MatrixOption::isElementMatrix<This>(), "[Error]: Sparse matrix should be element matrix");
    private:
        Utils::Array<ScalarType> elements;
        Utils::Array<size_t> minorIndexes;
        Utils::Array<size_t> majorStarts;
        size_t maxMinor;
    public:
        SparseMatrix();
        SparseMatrix(size_t row, size_t col);
        /* Operations */
        void insert(ScalarType x, size_t row, size_t col);
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const;
        [[nodiscard]] inline size_t getRow() const noexcept;
        [[nodiscard]] inline size_t getColumn() const noexcept;
        [[nodiscard]] inline size_t getMaxMajor() const noexcept;
        [[nodiscard]] inline size_t getMaxMinor() const noexcept;
        [[nodiscard]] size_t getNumNonZero() const noexcept { return elements.getLength(); }
    };

    template<class ScalarType, int option>
    SparseMatrix<ScalarType, option>::SparseMatrix()
            : Base()
            , elements()
            , minorIndexes()
            , majorStarts()
            , maxMinor(0) {}

    template<class ScalarType, int option>
    SparseMatrix<ScalarType, option>::SparseMatrix(size_t row, size_t col)
            : Base()
            , elements()
            , minorIndexes()
            , majorStarts((MatrixOption::isColumnMatrix<This>() ? col : row) + 1, 0)
            , maxMinor(MatrixOption::isColumnMatrix<This>() ? row : col) {}

    template<class ScalarType, int option>
    void SparseMatrix<ScalarType, option>::insert(ScalarType x, size_t row, size_t col) {
        assert(!x.isZero());
        assert(row < getRow());
        assert(col < getColumn());
        const size_t major = MatrixOption::isColumnMatrix<This>() ? col : row;
        const size_t minor = MatrixOption::isColumnMatrix<This>() ? row : col;
        /* Search existing element */ {
            const size_t from = majorStarts[major];
            const size_t to = majorStarts[major + 1];
            size_t index = 0;
            for (size_t i = from; i < to && index <= minor; ++i) {
                index = minorIndexes[i];
                if (index == minor) {
                    elements[i] = x;
                    return;
                }
            }
        }

        /* If not exist */ {
            for (size_t i = major + 1; i < majorStarts.getLength(); ++i)
                ++majorStarts[i];
            size_t insert_to = majorStarts[major];
            while (insert_to < minorIndexes.getLength() && minorIndexes[insert_to] < minor)
                ++insert_to;
            elements.insert(insert_to, x);
            minorIndexes.insert(insert_to, minor);
        }
    }

    template<class ScalarType, int option>
    ScalarType SparseMatrix<ScalarType, option>::calc(size_t row, size_t col) const {
        assert(row < getRow());
        assert(col < getColumn());
        const size_t major = MatrixOption::isColumnMatrix<This>() ? col : row;
        const size_t minor = MatrixOption::isColumnMatrix<This>() ? row : col;
        const size_t from = majorStarts[major];
        const size_t to = majorStarts[major + 1];
        size_t index = 0;
        for (size_t i = from; i < to && index <= minor; ++i) {
            index = minorIndexes[i];
            if (index == minor)
                return elements[i];
        }
        return ScalarType::Zero();
    }

    template<class ScalarType, int option>
    inline size_t SparseMatrix<ScalarType, option>::getRow() const noexcept {
        return MatrixOption::isColumnMatrix<This>() ? getMaxMinor() : getMaxMajor();
    }

    template<class ScalarType, int option>
    inline size_t SparseMatrix<ScalarType, option>::getColumn() const noexcept {
        return MatrixOption::isColumnMatrix<This>() ? getMaxMajor() : getMaxMinor();
    }

    template<class ScalarType, int option>
    inline size_t SparseMatrix<ScalarType, option>::getMaxMajor() const noexcept {
        return majorStarts.getLength() - 1;
    }

    template<class ScalarType, int option>
    inline size_t SparseMatrix<ScalarType, option>::getMaxMinor() const noexcept {
        return maxMinor;
    }
}
