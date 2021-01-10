/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

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
#ifndef PHYSICA_MATRIXIMPL_H
#define PHYSICA_MATRIXIMPL_H

namespace Physica::Core {
    /////////////////////////////////////ElementColumn Matrix//////////////////////////////////////////
    /*!
     * Warning: The elements have not been initialized, do not try to access them.
     */
    template<class T, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Element, maxRow, maxColumn>::DenseMatrix(size_t row, size_t column)
            : Utils::CStyleArray<T, maxRow * maxColumn>(row * column), row(row), column(column) {
        Q_ASSERT(row <= maxRow && column <= maxColumn);
    }
    /*!
     * Initializer list constructor for column matrix.
     *
     * e.g.
     * DenseMatrix<T, Column, maxRow, maxColumn> m{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}
     * We get a matrix looks like this:
     * (1, 4, 7)
     * (2, 5, 8)
     * (3, 6, 9)
     */
    template<class T, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Element, maxRow, maxColumn>::DenseMatrix(
            std::initializer_list<Utils::CStyleArray<T, maxRow>> list)
            : Utils::CStyleArray<T, maxRow * maxColumn>(list.begin()->getLength() * list.size())
            , row(list.begin()->getLength()), column(list.size()) {
        Q_ASSERT(row <= maxRow && column <= maxColumn);
        size_t i = 0;
        const auto end = list.end();
        for(auto ite = list.begin(); ite != end; ++ite, ++i) {
            const auto& array = *ite;
            for(size_t j = 0; j < row; ++j)
                Base::allocate(array[j], i * row + j);
        }
    }

    template<class T, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Element, maxRow, maxColumn>::DenseMatrix(
            DenseMatrix&& m) noexcept
            : Utils::CStyleArray<T, maxRow * maxColumn>(std::move(m)), row(m.row), column(m.column) {}

    template<class T, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Element, maxRow, maxColumn>&
    DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Element, maxRow, maxColumn>::operator=(
            DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Element, maxRow, maxColumn>&& m) noexcept {
        Base::operator=(static_cast<Base&&>(m));
        row = m.row;
        column = m.column;
    }
    /////////////////////////////////////ElementRow Matrix//////////////////////////////////////////
    /*!
    * Warning: The elements have not been initialized, do not try to access them.
    */
    template<class T, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Element, maxRow, maxColumn>::DenseMatrix(size_t row, size_t column)
            : Utils::CStyleArray<T, maxRow * maxColumn>(row * column), row(row), column(column) {
        Q_ASSERT(row <= maxRow && column <= maxColumn);
    }
    /*!
     * Initializer list constructor for row matrix.
     *
     * e.g.
     * DenseMatrix<T, Row, maxRow, maxColumn> m{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
     * We get a matrix looks like this:
     * (1, 2, 3)
     * (4, 5, 6)
     * (7, 8, 9)
     */
    template<class T, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Element, maxRow, maxColumn>::DenseMatrix(
            std::initializer_list<Utils::CStyleArray<T, maxColumn>> list)
                    : Utils::CStyleArray<T, maxRow * maxColumn>(list.begin()->getLength() * list.size())
                    , row(list.size()), column(list.begin()->getLength()) {
        Q_ASSERT(row <= maxRow && column <= maxColumn);
        size_t i = 0;
        const auto end = list.end();
        for(auto ite = list.begin(); ite != end; ++ite, ++i) {
            const auto& array = *ite;
            for(size_t j = 0; j < column; ++j)
                Base::allocate(array[j], i * column + j);
        }
    }

    template<class T, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Element, maxRow, maxColumn>::DenseMatrix(
            DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Element, maxRow, maxColumn>&& m) noexcept
            : Utils::CStyleArray<T, maxRow * maxColumn>(static_cast<Base&&>(m)), row(m.row), column(m.column) {}

    template<class T, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Element, maxRow, maxColumn>&
    DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Element, maxRow, maxColumn>::operator=(
            DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Element, maxRow, maxColumn>&& m) noexcept {
        Base::operator=(static_cast<Base&&>(m));
        row = m.row;
        column = m.column;
    }
    /////////////////////////////////////Utils::Dynamic Column Matrix//////////////////////////////////////////
    template<class T, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>::DenseMatrix(size_t length)
            : Utils::CStyleArray<Vector<T, maxRow>, maxColumn>(length) {}
            
    template<class T, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>::DenseMatrix(
            DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>&& m) noexcept
            : Utils::CStyleArray<Vector<T, maxRow>, maxColumn>(std::move(m)) {}

    template<class T, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>::DenseMatrix(
            std::initializer_list<VectorType> list)
            : Utils::CStyleArray<Vector<T, maxRow>, maxColumn>(list) {}

    template<class T, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>&
    DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>::operator=(
            DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>&& m) noexcept {
        Utils::CStyleArray<Vector<T, maxRow>, maxColumn>::operator=(std::move(m));
        return *this;
    }
    //!Optimize: may be inline append().
    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>::appendRow(const VectorType& v) {
        const auto c = getColumn();
        Q_ASSERT(c == v.getLength());
        for(size_t i = 0; i < c; ++i)
            Base::operator[](i).append(v[i]);
    }

    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>::appendRow(VectorType&& v) noexcept {
        const auto c = getColumn();
        Q_ASSERT(c == v.getLength());
        for(size_t i = 0; i < c; ++i)
            Base::operator[](i).append(std::move(v[i]));
    }

    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>::appendColumn(const VectorType& v) {
        Q_ASSERT(Base::getLength() == 0 || v.getLength() == getRow());
        append(v);
    }

    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>::appendColumn(VectorType&& v) noexcept {
        Q_ASSERT(Base::getLength() == 0 || v.getLength() == getRow());
        append(std::move(v));
    }

    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>::appendMatrixRow(const DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>& m) {
        const auto c = getColumn();
        Q_ASSERT(c == m.getColumn());
        for(size_t i = 0; i < c; ++i)
            Base::operator[](i).append(m[i]);
    }

    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>::appendMatrixRow(DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>&& m) {
        const auto c = getColumn();
        Q_ASSERT(c == m.getColumn());
        for(size_t i = 0; i < c; ++i)
            Base::operator[](i).append(std::move(m[i]));
    }

    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>::appendMatrixColumn(const DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>& m) {
        Q_ASSERT(Base::getLength() == 0 || getRow() == m.getRow());
        append(m);
    }

    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>::appendMatrixColumn(DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>&& m) {
        Q_ASSERT(Base::getLength() == 0 || getRow() == m.getRow());
        append(std::move(m));
    }

    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>::removeRowAt(size_t index) {
        for(size_t i = 0; i < getColumn(); ++i)
            Base::operator[](i).removeAt(index);
    }

    template<class T, size_t maxRow, size_t maxColumn>
    inline void DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>::removeColumnAt(size_t index) {
        Base::removeAt(index);
    }

    template<class T, size_t maxRow, size_t maxColumn>
    Vector<T, maxRow> DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>::cutRow() {
        const auto c = getColumn();
        VectorType result((Utils::CStyleArray<MultiScalar, Utils::Dynamic>(c)));
        for(size_t i = 0; i < c; ++i)
            result.allocate(Base::operator[](i).cutLast(), i);
        return result;
    }

    template<class T, size_t maxRow, size_t maxColumn>
    Vector<T, maxRow> DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>::cutColumn() {
        return Utils::CStyleArray<Vector<T, maxRow>, maxColumn>::cutLast();
    }

    template<class T, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, Utils::Dynamic, maxColumn>
            DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>::cutMatrixRow(size_t from) {
        const auto row = DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>::row();
        DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, Utils::Dynamic, maxColumn> result(row);
        for(size_t i = 0; i < row; ++i)
            result.allocate(VectorType(Base::operator[](i).cut(from)), i);
        return result;
    }

    template<class T, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, Utils::Dynamic>
            DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>::cutMatrixColumn(size_t from) {
        return DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, Utils::Dynamic>(Base::cut(from));
    }

    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>::rowSwap(size_t r1, size_t r2) noexcept {
        const auto length = getColumn();
        for(size_t i = 0; i < length; ++i) {
            auto& column = Base::operator[](i);
            Physica::Utils::swap(column[r1], column[r2]);
        }
    }

    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>::columnSwap(size_t c1, size_t c2) noexcept {
        Physica::Utils::swap(Base::operator[](c1), Base::operator[](c2));
    }
    //!Reduce the element at \r2 using \r1.
    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>::rowReduce(size_t r1, size_t r2, size_t element) {
        const auto& element_column = Base::operator[](element);
        Scalar dividend = element_column[r2] / element_column[r1];
        const auto length = getColumn();
        for(size_t i = 0; i < length; ++i) {
            auto& column = Base::operator[](i);
            column[r2] -= column[r1] * dividend;
        }
    }
    //!Reduce the element at \c2 using \c1.
    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>::columnReduce(size_t c1, size_t c2, size_t element) {
        Scalar dividend = operator()(element, c2) / operator()(element, c1);
        Base::operator[](c2) -= Base::operator[](c1) * dividend;
    }
    /////////////////////////////////////VectorRow Matrix//////////////////////////////////////////
    template<class T, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::DenseMatrix(size_t length)
            : Utils::CStyleArray<VectorType, Utils::Dynamic>(length) {}

    template<class T, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::DenseMatrix(std::initializer_list<VectorType> list)
            : Utils::CStyleArray<VectorType, Utils::Dynamic>(list) {}

    template<class T, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::DenseMatrix(
            DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>&& m) noexcept
            : Utils::CStyleArray<VectorType, Utils::Dynamic>(std::move(m)) {}

    template<class T, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>&
    DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::operator=(
            DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>&& m) noexcept {
        Utils::CStyleArray<VectorType, Utils::Dynamic>::operator=(std::move(m));
        return *this;
    }
    //!Optimize: may be inline append().
    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::appendRow(const VectorType& v) {
        Q_ASSERT(Base::getLength() == 0 || v.getLength() == getColumn());
        Base::append(v);
    }

    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::appendRow(VectorType&& v) noexcept {
        Q_ASSERT(Base::getLength() == 0 || v.getLength() == getColumn());
        Base::append(std::move(v));
    }

    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::appendColumn(const VectorType& v) {
        const auto length = (*this)[0].getLength();
        Q_ASSERT(length == v.getLength());
        for(size_t i = 0; i < length; ++i)
            Base::operator[](i).append(v[i]);
    }

    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::appendColumn(VectorType&& v) noexcept {
        const auto length = (*this)[0].getLength();
        Q_ASSERT(length == v.getLength());
        for(size_t i = 0; i < length; ++i)
            Base::operator[](i).append(std::move(v[i]));
    }

    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::appendMatrixRow(const DenseMatrix& m) {
        Q_ASSERT(Base::operator[](0).getLength() == m[0].getLength());
        append(m);
    }

    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::appendMatrixRow(DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>&& m) {
        Q_ASSERT(Base::operator[](0).getLength() == m[0].getLength());
        append(std::move(m));
    }

    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::appendMatrixColumn(DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>&& m) {
        const auto length = getRow();
        Q_ASSERT(length == m.getLength());
        for(size_t i = 0; i < length; ++i)
            Base::operator[](i).append(std::move(m[i]));
    }

    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::appendMatrixColumn(const DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>& m) {
        const auto length = getRow();
        Q_ASSERT(length == m.getLength());
        for(size_t i = 0; i < length; ++i)
            Base::operator[](i).append(m[i]);
    }

    template<class T, size_t maxRow, size_t maxColumn>
    inline void DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::removeRowAt(size_t index) {
        Base::removeAt(index);
    }

    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::removeColumnAt(size_t index) {
        for(size_t i = 0; i < getRow(); ++i)
            Base::operator[](i).removeAt(index);
    }

    template<class T, size_t maxRow, size_t maxColumn>
    Vector<T, maxColumn> DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::cutRow() {
        return Base::cutLast();
    }

    template<class T, size_t maxRow, size_t maxColumn>
    Vector<T, maxColumn> DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::cutColumn() {
        const auto row = DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::row();
        VectorType result((Utils::CStyleArray<MultiScalar, Utils::Dynamic>(row)));
        for(size_t i = 0; i < row; ++i)
            result.allocate((*this)[i].cutLast(), i);
        return result;
    }

    template<class T, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, Utils::Dynamic, maxColumn>
            DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::cutMatrixRow(size_t from) {
        return DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, Utils::Dynamic, maxColumn>(Base::cut(from));
    }

    template<class T, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, Utils::Dynamic>
            DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::cutMatrixColumn(size_t from) {
        const auto column = DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::column();
        DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, Utils::Dynamic> result(column);
        for(size_t i = 0; i < column; ++i)
            result.allocate(VectorType((*this)[i].cut(from)), i);
        return result;
    }

    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::rowSwap(size_t r1, size_t r2) noexcept {
        Physica::Utils::swap((*this)[r1], (*this)[r2]);
    }

    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::columnSwap(size_t c1, size_t c2) noexcept {
        const auto length = getRow();
        for(size_t i = 0; i < length; ++i) {
            auto& row = (*this)[i];
            Physica::Utils::swap(row[c1], row[c2]);
        }
    }
    //!Reduce the element at \r2 using \r1
    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::rowReduce(size_t r1, size_t r2, size_t element) {
        Scalar dividend = (*this)(element, r2) / (*this)(element, r1);
        (*this)[r2] -= (*this)[r1] * dividend;
    }
    //!Reduce the element at \c2 using \c1.
    template<class T, size_t maxRow, size_t maxColumn>
    void DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>::columnReduce(size_t c1, size_t c2, size_t element) {
        const auto& element_row = Base::operator[](element);
        Scalar dividend = element_row[c2] / element_row[c1];
        const auto length = getRow();
        for(size_t i = 0; i < length; ++i) {
            auto& row = (*this)[i];
            row[c2] -= row[c1] * dividend;
        }
    }
    ////////////////////////////////////////Global////////////////////////////////////////////
    /*!
     * Print all elements.
     */
    template<class T, int type, size_t maxRow, size_t maxColumn>
    std::ostream& operator<<(std::ostream& os, const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto row = m.getRow();
        const auto column = m.getColumn();
        //10 is the max precision of double.
        os << std::setprecision(10);
        for(size_t i = 0; i < row; ++i) {
            for(size_t j = 0; j < column; ++j)
                os << m(i, j) << '\t';
            os << '\n';
        }
        //6 is the default precision.
        return os << std::setprecision(6);
    }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> operator+(
            const DenseMatrix<T, type, maxRow, maxColumn>& m1, const DenseMatrix<T, type, maxRow, maxColumn>& m2) {
        Q_ASSERT(m1.getRow() == m2.getRow() && m1.getColumn() == m2.getColumn());
        const auto length = m1.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(m1[i] + m2[i], i);
        return result;
    }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> operator-(
        const DenseMatrix<T, type, maxRow, maxColumn>& m1, const DenseMatrix<T, type, maxRow, maxColumn>& m2) {
        Q_ASSERT(m1.getRow() == m2.getRow() && m1.getColumn() == m2.getColumn());
        const auto length = m1.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(m1[i] - m2[i], i);
        return result;
    }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> operator*(
            const DenseMatrix<T, type, maxRow, maxColumn>& m1, const DenseMatrix<T, type, maxRow, maxColumn>& m2) {
        Q_ASSERT(m1.getColumn() == m2.getRow());
        constexpr bool isColumn = type == DenseMatrixType::Column | DenseMatrixType::Element || type == DenseMatrixType::Column | DenseMatrixType::Vector;

        const auto result_row = m1.row();
        const auto result_column = m2.column();
        const auto m1_column = m1.column();
        const auto matrix_size = isColumn ? result_column : result_row;
        const auto vector_length = isColumn ? result_row : result_column;
        DenseMatrix<T, type, maxRow, maxColumn> result(matrix_size);

        size_t i, j;
        const size_t& r = isColumn ? j : i;
        const size_t& c = isColumn ? i : j;
        for(i = 0; i < matrix_size; ++i) {
            Vector<T, isColumn ? maxRow : maxColumn> new_vector(vector_length);
            for(j = 0; j < vector_length; ++j) {
                T element = T::getZero();
                for(size_t k = 0; k < m1_column; ++k)
                    element += m1(r, k) * m2(k, c);
                new_vector.allocate(std::move(element), j);
            }
        }
        return result;
    }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> operator*(
        const DenseMatrix<T, type, maxRow, maxColumn>& m, const MultiScalar& n) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(m[i] * n, i);
        return result;
    }
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> reciprocal(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(reciprocal(m[i]), i);
        return result;
    }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> sqrt(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(sqrt(m[i]), i);
        return result;
    }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> factorial(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(factorial(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> ln(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(ln(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> log(const DenseMatrix<T, type, maxRow, maxColumn>& m, const MultiScalar& a) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(log(m[i], a), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> exp(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(exp(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> pow(const DenseMatrix<T, type, maxRow, maxColumn>& m, const MultiScalar& a) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(pow(m[i], a), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> cos(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(cos(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> sin(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(sin(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> tan(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(tan(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> sec(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(sec(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> csc(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(csc(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> cot(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(cot(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arccos(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arccos(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arcsin(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arcsin(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arctan(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arctan(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arcsec(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arcsec(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arccsc(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arccsc(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arccot(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arccot(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> cosh(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(cosh(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> sinh(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(sinh(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> tanh(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(tanh(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> sech(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(sech(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> csch(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(csch(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> coth(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(coth(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arccosh(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arccosh(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arcsinh(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arcsinh(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arctanh(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arctanh(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arcsech(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arcsech(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arccsch(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arccsch(m[i]), i);
        return result;
    }
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arccoth(const DenseMatrix<T, type, maxRow, maxColumn>& m) {
        const auto length = m.getLength();
        DenseMatrix<T, type, maxRow, maxColumn> result(length);
        for(size_t i = 0; i < length; ++i)
            result.allocate(arccoth(m[i]), i);
        return result;
    }
}

#endif