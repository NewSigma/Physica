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

namespace Physica::Core::Internal {
    //////////////////////////////////////////////Column-Element//////////////////////////////////////////////
    template<class Derived>
    template<size_t Length, size_t MaxLength>
    void AbstractDenseMatrixStorage<Derived, DenseMatrixType::Column | DenseMatrixType::Element>::appendRow(
            const Vector<T, Length, MaxLength>& v) {
        constexpr size_t MaxRowAtCompile = Traits<Derived>::MaxRowAtCompile;
        static_assert(Traits<Derived>::RowAtCompile == Dynamic);
        assert(MaxRowAtCompile == Dynamic || MaxRowAtCompile > getDerived().getRow());
        assert(v.getLength() == getDerived().getColumn());

        const size_t column = getDerived().getColumn();
        const size_t newSize = Base::getLength() + column;
        if (newSize > Base::getCapacity())
            Base::reserve(newSize);
        
        const size_t row = getDerived().getRow();
        /* Move elements */ {
            T* data = Base::data();
            const size_t row_1 = row + 1;
            for (size_t i = column - 1; i > 1; --i)
                memmove(row_1 * i, row * i, row * sizeof(T));
        }
        /* Insert elements */ {
            size_t index = row; //Insert first element to position row.
            const size_t newRow = row + 1;
            for (auto ite = v.begin(); ite != v.end(); ++ite) {
                Base::init(*ite, index);
                index += newRow;
            }
            Base::setLength(newSize);
        }
    }

    template<class Derived>
    void AbstractDenseMatrixStorage<Derived, DenseMatrixType::Column | DenseMatrixType::Element>::removeColumnAt(size_t index) {
        Derived& matrix = getDerived();
        assert(index < matrix.getColumn());
    }

    template<class Derived>
    void AbstractDenseMatrixStorage<Derived, DenseMatrixType::Column | DenseMatrixType::Element>::rowSwap(size_t r1, size_t r2) {
        Derived& matrix = getDerived();
        const size_t row = matrix.getRow();
        const size_t column = matrix.getColumn();
        assert(r1 < row && r2 < row);
        for (size_t i = 0, temp = 0; i < column; ++i, temp += column)
            swap(matrix[temp + r1], matrix[temp + r2]);
    }

    template<class Derived>
    void AbstractDenseMatrixStorage<Derived, DenseMatrixType::Column | DenseMatrixType::Element>::columnSwap(size_t c1, size_t c2) {
        Derived& matrix = getDerived();
        const size_t row = matrix.getRow();
        const size_t column = matrix.getColumn();
        assert(c1 < column && c2 < column);
        const size_t offset1 = c1 * row;
        const size_t offset2 = c2 * row;
        for (size_t i = 0; i < row; ++i) {
            swap(matrix[offset1 + i], matrix[offset2 + i]);
        }
    }
    //////////////////////////////////////////////Row-Element//////////////////////////////////////////////
    template<class Derived>
    template<size_t Length, size_t MaxLength>
    void AbstractDenseMatrixStorage<Derived, DenseMatrixType::Row | DenseMatrixType::Element>::appendRow(
            const Vector<T, Length, MaxLength>& v) {
        constexpr size_t MaxRowAtCompile = Traits<Derived>::MaxRowAtCompile;
        static_assert(Traits<Derived>::RowAtCompile == Dynamic);
        assert(MaxRowAtCompile == Dynamic || MaxRowAtCompile > getDerived().getRow());
        assert(v.getLength() == getDerived().getColumn());

        const size_t column = getDerived().getColumn();
        const size_t length = Base::getLength();
        const size_t newSize = length + column;
        if (newSize > Base::getCapacity())
            Base::reserve(newSize);
        /* Insert elements */ {
            for (size_t i = length, j = 0; i < newSize; ++i, ++j)
                Base::init(v[j], j);
            Base::setLength(newSize);
        }
    }

    template<class Derived>
    void AbstractDenseMatrixStorage<Derived, DenseMatrixType::Row | DenseMatrixType::Element>::removeColumnAt(size_t index) {
        Derived& matrix = getDerived();
        assert(index < matrix.getColumn());
    }

    template<class Derived>
    void AbstractDenseMatrixStorage<Derived, DenseMatrixType::Row | DenseMatrixType::Element>::rowSwap(size_t r1, size_t r2) {
        Derived& matrix = getDerived();
        const size_t row = matrix.getRow();
        const size_t column = matrix.getColumn();
        assert(r1 < row && r2 < row);
        const size_t offset1 = r1 * column;
        const size_t offset2 = r2 * column;
        for (size_t i = 0; i < column; ++i)
            swap(matrix[offset1 + i], matrix[offset2 + i]);
    }

    template<class Derived>
    void AbstractDenseMatrixStorage<Derived, DenseMatrixType::Row | DenseMatrixType::Element>::columnSwap(size_t c1, size_t c2) {
        Derived& matrix = getDerived();
        const size_t row = matrix.getRow();
        const size_t column = matrix.getColumn();
        assert(c1 < column && c2 < column);
        for (size_t i = 0, temp = 0; i < column; ++i, temp += row)
            swap(matrix[temp + c1], matrix[temp + c2]);
    }
    //////////////////////////////////////////////Column-Vector//////////////////////////////////////////////
    template<class Derived>
    template<size_t Length, size_t MaxLength>
    void AbstractDenseMatrixStorage<Derived, DenseMatrixType::Column | DenseMatrixType::Vector>::appendRow(
            const Vector<T, Length, MaxLength>& v) {
        constexpr size_t MaxRowAtCompile = Traits<Derived>::MaxRowAtCompile;
        static_assert(Traits<Derived>::RowAtCompile == Dynamic);
        assert(MaxRowAtCompile == Dynamic || MaxRowAtCompile > getDerived().getRow());
        assert(v.getLength() == getDerived().getColumn());

        const auto end = Base::end();
        auto ite = Base::begin();
        auto ite1 = v.begin();
        for (; ite != end(); ++ite, ++ite1)
            (*ite).append(*ite1);
    }

    template<class Derived>
    void AbstractDenseMatrixStorage<Derived, DenseMatrixType::Column | DenseMatrixType::Vector>::removeColumnAt(size_t index) {
        Derived& matrix = getDerived();
        assert(index < matrix.getColumn());
    }

    template<class Derived>
    void AbstractDenseMatrixStorage<Derived, DenseMatrixType::Column | DenseMatrixType::Vector>::rowSwap(size_t r1, size_t r2) {
        Derived& matrix = getDerived();
        const size_t row = matrix.getRow();
        assert(r1 < row && r2 < row);
        for (auto& columnVector : matrix)
            columnVector[r1].swap(columnVector[r2]);
    }

    template<class Derived>
    void AbstractDenseMatrixStorage<Derived, DenseMatrixType::Column | DenseMatrixType::Vector>::columnSwap(size_t c1, size_t c2) {
        Derived& matrix = getDerived();
        [[maybe_unused]] const size_t column = matrix.getColumn();
        assert(c1 < column && c2 < column);
        matrix[c1].swap(matrix[c2]);
    }
    //////////////////////////////////////////////Row-Vector//////////////////////////////////////////////
    template<class Derived>
    template<size_t Length, size_t MaxLength>
    void AbstractDenseMatrixStorage<Derived, DenseMatrixType::Row | DenseMatrixType::Vector>::appendRow(
            const Vector<T, Length, MaxLength>& v) {
        constexpr size_t MaxRowAtCompile = Traits<Derived>::MaxRowAtCompile;
        static_assert(Traits<Derived>::RowAtCompile == Dynamic);
        assert(MaxRowAtCompile == Dynamic || MaxRowAtCompile > getDerived().getRow());
        assert(v.getLength() == getDerived().getColumn());

        Base::append(v);
    }

    template<class Derived>
    void AbstractDenseMatrixStorage<Derived, DenseMatrixType::Row | DenseMatrixType::Vector>::removeColumnAt(size_t index) {
        Derived& matrix = getDerived();
        assert(index < matrix.getColumn());
    }

    template<class Derived>
    void AbstractDenseMatrixStorage<Derived, DenseMatrixType::Row | DenseMatrixType::Vector>::rowSwap(size_t r1, size_t r2) {
        Derived& matrix = getDerived();
        [[maybe_unused]] const size_t row = matrix.getRow();
        assert(r1 < row && r2 < row);
        matrix[r1].swap(matrix[r2]);
    }

    template<class Derived>
    void AbstractDenseMatrixStorage<Derived, DenseMatrixType::Row | DenseMatrixType::Vector>::columnSwap(size_t c1, size_t c2) {
        Derived& matrix = getDerived();
        const size_t column = matrix.getColumn();
        assert(c1 < column && c2 < column);
        for (auto& rowVector : matrix)
            rowVector[c1].swap(rowVector[c2]);
    }
}
