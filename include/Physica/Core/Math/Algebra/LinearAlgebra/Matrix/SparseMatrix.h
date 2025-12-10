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

namespace Physica {
    /**
     * CSR & CSC format
     */
    template<Scalar T, int Option = MatrixOption::Row>
    class SparseMatrix : public RValueMatrix<SparseMatrix<T, Option>> {
        using This = SparseMatrix<T, Option>;
        using Base = RValueMatrix<This>;
    private:
        Array<T> elements;
        Array<size_t> minorIndexes;
        Array<size_t> majorStarts;
        size_t maxMinor = 0;
    public:
        SparseMatrix() = default;
        SparseMatrix(size_t row, size_t col);
        SparseMatrix(const This&) = default;
        SparseMatrix(This&&) noexcept = default;
        ~SparseMatrix() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void insert(T x, size_t row, size_t col);
        void resize(size_t row, size_t col);
        void clear();
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] T calc(size_t row, size_t col) const;
        [[nodiscard]] size_t getRow() const noexcept;
        [[nodiscard]] size_t getCol() const noexcept;
        [[nodiscard]] size_t getMaxMajor() const noexcept;
        [[nodiscard]] size_t getMaxMinor() const noexcept;
        [[nodiscard]] size_t getNumNonZero() const noexcept { return elements.getLength(); }
        [[nodiscard]] const auto& getElements() const { return elements; }
        [[nodiscard]] const auto& getMinorIndexes() const { return minorIndexes; }
        [[nodiscard]] const auto& getMajorStarts() const { return majorStarts; }
    };

    template<Scalar T, int Option>
    SparseMatrix<T, Option>::SparseMatrix(size_t row, size_t col)
            : Base()
            , elements()
            , majorStarts(MatrixOption::selectMajor<This>(row, col) + 1, 0)
            , maxMinor(MatrixOption::selectMinor<This>(row, col)) {
        size_t size = std::max(row, col);
        elements.reserve(size);
        minorIndexes.reserve(size);
    }

    template<Scalar T, int Option>
    void SparseMatrix<T, Option>::insert(T x, size_t row, size_t col) {
        // TODO: Inserting 0 element is useless
        assert(row < getRow());
        assert(col < getCol());
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

    template<Scalar T, int Option>
    void SparseMatrix<T, Option>::resize(size_t row, size_t col) {
        elements.resize(0);
        minorIndexes.resize(0);
        majorStarts.resize(MatrixOption::selectMajor<This>(row, col) + 1, 0);
        maxMinor = MatrixOption::selectMinor<This>(row, col);
    }

    template<Scalar T, int Option>
    void SparseMatrix<T, Option>::clear() {
        elements.resize(0);
        minorIndexes.resize(0);
        for (auto& i : majorStarts)
            i = 0;
    }

    template<Scalar T, int Option>
    void SparseMatrix<T, Option>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        elements.swap(obj.elements);
        minorIndexes.swap(obj.minorIndexes);
        majorStarts.swap(obj.majorStarts);
        std::swap(maxMinor, obj.maxMinor);
    }

    template<Scalar T, int Option>
    T SparseMatrix<T, Option>::calc(size_t row, size_t col) const {
        assert(row < getRow());
        assert(col < getCol());
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
        return T(0);
    }

    template<Scalar T, int Option>
    size_t SparseMatrix<T, Option>::getRow() const noexcept {
        return MatrixOption::isColMatrix<This>() ? getMaxMinor() : getMaxMajor();
    }

    template<Scalar T, int Option>
    size_t SparseMatrix<T, Option>::getCol() const noexcept {
        return MatrixOption::isColMatrix<This>() ? getMaxMajor() : getMaxMinor();
    }

    template<Scalar T, int Option>
    size_t SparseMatrix<T, Option>::getMaxMajor() const noexcept {
        return majorStarts.getLength() - 1;
    }

    template<Scalar T, int Option>
    size_t SparseMatrix<T, Option>::getMaxMinor() const noexcept {
        return maxMinor;
    }
}

namespace Physica {
    template<Scalar T, int Op>
    class Traits<SparseMatrix<T, Op>> {
    public:
        using ScalarType = T;
        constexpr static int Option = Op;
        constexpr static size_t RowAtCompile = Dynamic;
        constexpr static size_t ColAtCompile = Dynamic;
        constexpr static size_t SizeAtCompile = Dynamic;
    };
}

#include "MatrixImpl/MatrixProduct/SpMV.h"
