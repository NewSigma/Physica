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
        for (size_t i = 0; i < column; ++i) {
            size_t temp = i * column;
            swap(temp + r1, temp + r2);
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
            swap(offset1 + i, offset2 + i);
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
            swap(columnVector[r1], columnVector[r2]);
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
        using std::swap;
        Derived& matrix = getDerived();
        const size_t row = matrix.getRow();
        assert(r1 < row && r2 < row);
        swap(matrix[r1], matrix[r2]);
    }
}
