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

#include "MatrixImpl/RValueMatrix.h"
#include "MatrixOption.h"

namespace Physica::Core {
    template<class ScalarType, int Option = MatrixOption::Row | MatrixOption::Element>
    class SparseMatrix : public RValueMatrix<SparseMatrix<ScalarType, Option>> {
        using This = SparseMatrix<ScalarType, Option>;
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
        void resize(size_t row, size_t col);
        void clear();
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const;
        [[nodiscard]] inline size_t getRow() const noexcept;
        [[nodiscard]] inline size_t getColumn() const noexcept;
        [[nodiscard]] const Utils::Array<ScalarType>& getElements() const { return elements; }
        [[nodiscard]] const Utils::Array<size_t>& getMinorIndexes() const { return minorIndexes; }
        [[nodiscard]] const Utils::Array<size_t>& getMajorStarts() const { return majorStarts; }
        [[nodiscard]] inline size_t getMaxMajor() const noexcept;
        [[nodiscard]] inline size_t getMaxMinor() const noexcept;
        [[nodiscard]] size_t getNumNonZero() const noexcept { return elements.getLength(); }
    };

    template<class ScalarType, int Option>
    SparseMatrix<ScalarType, Option>::SparseMatrix()
            : Base()
            , elements()
            , minorIndexes()
            , majorStarts()
            , maxMinor(0) {}

    template<class ScalarType, int Option>
    SparseMatrix<ScalarType, Option>::SparseMatrix(size_t row, size_t col)
            : Base()
            , elements()
            , minorIndexes()
            , majorStarts(MatrixOption::selectMajor<This>(row, col) + 1, 0)
            , maxMinor(MatrixOption::selectMinor<This>(row, col)) {}

    template<class ScalarType, int Option>
    void SparseMatrix<ScalarType, Option>::insert(ScalarType x, size_t row, size_t col) {
        // TODO: Inserting 0 element is useless
        assert(row < getRow());
        assert(col < getColumn());
        const size_t major = MatrixOption::selectMajor<This>(row, col);
        const size_t minor = MatrixOption::selectMinor<This>(row, col);
        const size_t from = majorStarts[major];
        const size_t to = majorStarts[major + 1];
        /* Search existing element */ {
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
            size_t insert_to = from;
            while (insert_to < to && minorIndexes[insert_to] < minor)
                ++insert_to;
            elements.insert(insert_to, x);
            minorIndexes.insert(insert_to, minor);
            for (size_t i = major + 1; i < majorStarts.getLength(); ++i)
                ++majorStarts[i];
        }
    }

    template<class ScalarType, int Option>
    void SparseMatrix<ScalarType, Option>::resize(size_t row, size_t col) {
        elements.resize(0);
        minorIndexes.resize(0);
        majorStarts.resize(MatrixOption::selectMajor<This>(row, col) + 1, 0);
        maxMinor = MatrixOption::selectMinor<This>(row, col);
    }

    template<class ScalarType, int Option>
    void SparseMatrix<ScalarType, Option>::clear() {
        elements.resize(0);
        minorIndexes.resize(0);
        for (auto& i : majorStarts)
            i = 0;
    }

    template<class ScalarType, int Option>
    ScalarType SparseMatrix<ScalarType, Option>::calc(size_t row, size_t col) const {
        assert(row < getRow());
        assert(col < getColumn());
        const size_t major = MatrixOption::selectMajor<This>(row, col);
        const size_t minor = MatrixOption::selectMinor<This>(row, col);
        const size_t from = majorStarts[major];
        const size_t to = majorStarts[major + 1];
        size_t index = 0;
        for (size_t i = from; i < to && index <= minor; ++i) {
            index = minorIndexes[i];
            if (index == minor)
                return elements[i];
        }
        return ScalarType(0);
    }

    template<class ScalarType, int Option>
    inline size_t SparseMatrix<ScalarType, Option>::getRow() const noexcept {
        return MatrixOption::isColumnMatrix<This>() ? getMaxMinor() : getMaxMajor();
    }

    template<class ScalarType, int Option>
    inline size_t SparseMatrix<ScalarType, Option>::getColumn() const noexcept {
        return MatrixOption::isColumnMatrix<This>() ? getMaxMajor() : getMaxMinor();
    }

    template<class ScalarType, int Option>
    inline size_t SparseMatrix<ScalarType, Option>::getMaxMajor() const noexcept {
        return majorStarts.getLength() - 1;
    }

    template<class ScalarType, int Option>
    inline size_t SparseMatrix<ScalarType, Option>::getMaxMinor() const noexcept {
        return maxMinor;
    }
}

namespace Physica {
    template<class T, int Op>
    class Traits<Core::SparseMatrix<T, Op>> {
    public:
        using ScalarType = T;
        constexpr static int Option = Op;
        constexpr static size_t RowAtCompile = Dynamic;
        constexpr static size_t ColumnAtCompile = Dynamic;
        constexpr static size_t MaxRowAtCompile = Dynamic;
        constexpr static size_t MaxColumnAtCompile = Dynamic;
        constexpr static size_t SizeAtCompile = Dynamic;
        constexpr static size_t MaxSizeAtCompile = Dynamic;
    };
}

#include "MatrixProduct/SparseMatrixProduct.h"
