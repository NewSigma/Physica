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

#include "AbstractDenseMatrixStorage.h"

namespace Physica::Core {
    using Utils::Dynamic;
    /**
     * This class handles specifications of \tparam Row, \param Column, \tparam MaxRow and \tparam MaxColumn.
     * 
     * Optimize: May be DenseMatrixStorage take only one
     * template argument class T and move implementation to Internal?
     */
    template<class T, int type, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage;

    namespace Internal {
        using namespace Utils;

        template<class T>
        class Traits;

        template<class T, int type, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
        class Traits<DenseMatrixStorage<T, type, Row, Column, MaxRow, MaxColumn>> {
        public:
            using ScalarType = T;
            constexpr static int MatrixOption = type;
            constexpr static size_t RowAtCompile = Row;
            constexpr static size_t ColumnAtCompile = Column;
            constexpr static size_t MaxRowAtCompile = MaxRow;
            constexpr static size_t MaxColumnAtCompile = MaxColumn;
            constexpr static size_t SizeAtCompile = Row * Column;
            constexpr static size_t MaxSizeAtCompile = MaxRow * MaxColumn;
        };
    }
    ////////////////////////////////////ColumnElement////////////////////////////////////
    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixOption::Column | DenseMatrixOption::Element, Row, Column, MaxRow, MaxColumn>
            : public Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Column | DenseMatrixOption::Element
                                                                           , Row
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Column | DenseMatrixOption::Element> {
        using Base = Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Column | DenseMatrixOption::Element
                                                                           , Row
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Column | DenseMatrixOption::Element>;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return Row; }
        [[nodiscard]] constexpr static size_t getColumn() noexcept { return Column; }
    };

    template<class T, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixOption::Column | DenseMatrixOption::Element, Dynamic, Column, MaxRow, MaxColumn>
            : public Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Column | DenseMatrixOption::Element
                                                                           , Dynamic
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Column | DenseMatrixOption::Element> {
        using Base = Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Column | DenseMatrixOption::Element
                                                                           , Dynamic
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Column | DenseMatrixOption::Element>;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { assert(Base::getLength() % getColumn() == 0); return Base::getLength() / getColumn(); }
        [[nodiscard]] constexpr static size_t getColumn() noexcept { return Column; }
    };

    template<class T, size_t Row, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixOption::Column | DenseMatrixOption::Element, Row, Dynamic, MaxRow, MaxColumn>
            : public Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Column | DenseMatrixOption::Element
                                                                           , Row
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Column | DenseMatrixOption::Element> {
        using Base = Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Column | DenseMatrixOption::Element
                                                                           , Row
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Column | DenseMatrixOption::Element>;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return Row; }
        [[nodiscard]] size_t getColumn() const noexcept { assert(Base::getLength() % getRow() == 0); return Base::getLength() / getRow(); }
    };

    template<class T, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixOption::Column | DenseMatrixOption::Element, Dynamic, Dynamic, MaxRow, MaxColumn>
            : public Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Column | DenseMatrixOption::Element
                                                                           , Dynamic
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Column | DenseMatrixOption::Element> {
        using Base = Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Column | DenseMatrixOption::Element
                                                                           , Dynamic
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Column | DenseMatrixOption::Element>;
    private:
        size_t row;
    public:
        DenseMatrixStorage() : Base(), row(0) {}
        DenseMatrixStorage(size_t row_, size_t column) : Base(row_, column), row(row_) {}
        DenseMatrixStorage(size_t row_, size_t column, const T& t) : Base(row_, column, t), row(row_) {}
        DenseMatrixStorage(size_t row_, std::initializer_list<T> list) : Base(list), row(row_) { assert(list.size() % row == 0); }
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return row; }
        [[nodiscard]] size_t getColumn() const noexcept { assert(Base::getLength() % getRow() == 0); return Base::getLength() / getRow(); }
    protected:
        /* Helper */
        void swap(DenseMatrixStorage& storage) noexcept { Base::swap(storage); std::swap(row, storage.row); }
    };
    ////////////////////////////////////RowElement////////////////////////////////////
    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixOption::Row | DenseMatrixOption::Element, Row, Column, MaxRow, MaxColumn>
            : public Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Row | DenseMatrixOption::Element
                                                                           , Row
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Row | DenseMatrixOption::Element> {
        using Base = Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Row | DenseMatrixOption::Element
                                                                           , Row
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Row | DenseMatrixOption::Element>;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return Row; }
        [[nodiscard]] constexpr static size_t getColumn() noexcept { return Column; }
    };

    template<class T, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixOption::Row | DenseMatrixOption::Element, Dynamic, Column, MaxRow, MaxColumn>
            : public Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Row | DenseMatrixOption::Element
                                                                           , Dynamic
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Row | DenseMatrixOption::Element> {
        using Base = Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Row | DenseMatrixOption::Element
                                                                           , Dynamic
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Row | DenseMatrixOption::Element>;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { assert(Base::getLength() % getColumn() == 0); return Base::getLength() / getColumn(); }
        [[nodiscard]] constexpr static size_t getColumn() noexcept { return Column; }
    };

    template<class T, size_t Row, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixOption::Row | DenseMatrixOption::Element, Row, Dynamic, MaxRow, MaxColumn>
            : public Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Row | DenseMatrixOption::Element
                                                                           , Row
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Row | DenseMatrixOption::Element> {
        using Base = Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Row | DenseMatrixOption::Element
                                                                           , Row
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Row | DenseMatrixOption::Element>;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return Row; }
        [[nodiscard]] size_t getColumn() const noexcept { assert(Base::getLength() % getRow() == 0); return Base::getLength() / getRow(); }
    };

    template<class T, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixOption::Row | DenseMatrixOption::Element, Dynamic, Dynamic, MaxRow, MaxColumn>
            : public Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Row | DenseMatrixOption::Element
                                                                           , Dynamic
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Row | DenseMatrixOption::Element> {
        using Base = Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Row | DenseMatrixOption::Element
                                                                           , Dynamic
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Row | DenseMatrixOption::Element>;
    private:
        size_t row;
    public:
        DenseMatrixStorage() : Base(), row(0) {}
        DenseMatrixStorage(size_t row_, size_t column) : Base(row_, column), row(row_) {}
        DenseMatrixStorage(size_t row_, size_t column, const T& t) : Base(row_, column, t), row(row_) {}
        DenseMatrixStorage(size_t row_, std::initializer_list<T> list) : Base(list), row(row_) { assert(list.size() % row == 0); }
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return row; }
        [[nodiscard]] size_t getColumn() const noexcept { assert(Base::getLength() % getRow() == 0); return Base::getLength() / getRow(); }
    protected:
        /* Helper */
        void swap(DenseMatrixStorage& storage) noexcept { Base::swap(storage); std::swap(row, storage.row); }
    };
    ////////////////////////////////////ColumnVector////////////////////////////////////
    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixOption::Column | DenseMatrixOption::Vector, Row, Column, MaxRow, MaxColumn>
            : public Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Column | DenseMatrixOption::Vector
                                                                           , Row
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Column | DenseMatrixOption::Vector> {
        using Base = Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Column | DenseMatrixOption::Vector
                                                                           , Row
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Column | DenseMatrixOption::Vector>;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return Row; }
        [[nodiscard]] constexpr static size_t getColumn() noexcept { return Column; }
    };

    template<class T, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixOption::Column | DenseMatrixOption::Vector, Dynamic, Column, MaxRow, MaxColumn>
            : public Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Column | DenseMatrixOption::Vector
                                                                           , Dynamic
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Column | DenseMatrixOption::Vector> {
        using Base = Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Column | DenseMatrixOption::Vector
                                                                           , Dynamic
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Column | DenseMatrixOption::Vector>;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return Base::operator[](0).getLength(); }
        [[nodiscard]] constexpr static size_t getColumn() noexcept { return Column; }
    };

    template<class T, size_t Row, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixOption::Column | DenseMatrixOption::Vector, Row, Dynamic, MaxRow, MaxColumn>
            : public Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Column | DenseMatrixOption::Vector
                                                                           , Row
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Column | DenseMatrixOption::Vector> {
        using Base = Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Column | DenseMatrixOption::Vector
                                                                           , Row
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Column | DenseMatrixOption::Vector>;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return Row; }
        [[nodiscard]] size_t getColumn() const noexcept {return Base::getLength(); }
    };

    template<class T, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixOption::Column | DenseMatrixOption::Vector, Dynamic, Dynamic, MaxRow, MaxColumn>
            : public Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Column | DenseMatrixOption::Vector
                                                                           , Dynamic
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Column | DenseMatrixOption::Vector> {
        using Base = Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Column | DenseMatrixOption::Vector
                                                                           , Dynamic
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Column | DenseMatrixOption::Vector>;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return Base::operator[](0).getLength(); }
        [[nodiscard]] size_t getColumn() const noexcept { return Base::getLength(); }
    };
    ////////////////////////////////////RowVector////////////////////////////////////
    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixOption::Row | DenseMatrixOption::Vector, Row, Column, MaxRow, MaxColumn>
            : public Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Row | DenseMatrixOption::Vector
                                                                           , Row
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Row | DenseMatrixOption::Vector> {
        using Base = Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Row | DenseMatrixOption::Vector
                                                                           , Row
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Row | DenseMatrixOption::Vector>;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return Row; }
        [[nodiscard]] constexpr static size_t getColumn() noexcept { return Column; }
    };

    template<class T, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixOption::Row | DenseMatrixOption::Vector, Dynamic, Column, MaxRow, MaxColumn>
            : public Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Row | DenseMatrixOption::Vector
                                                                           , Dynamic
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Row | DenseMatrixOption::Vector> {
        using Base = Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Row | DenseMatrixOption::Vector
                                                                           , Dynamic
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Row | DenseMatrixOption::Vector>;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return Base::getLength(); }
        [[nodiscard]] constexpr static size_t getColumn() noexcept { return Column; }
    };

    template<class T, size_t Row, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixOption::Row | DenseMatrixOption::Vector, Row, Dynamic, MaxRow, MaxColumn>
            : public Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Row | DenseMatrixOption::Vector
                                                                           , Row
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Row | DenseMatrixOption::Vector> {
        using Base = Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Row | DenseMatrixOption::Vector
                                                                           , Row
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Row | DenseMatrixOption::Vector>;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return Row; }
        [[nodiscard]] size_t getColumn() const noexcept {return Base::operator[](0).getLength(); }
    };

    template<class T, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixOption::Row | DenseMatrixOption::Vector, Dynamic, Dynamic, MaxRow, MaxColumn>
            : public Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Row | DenseMatrixOption::Vector
                                                                           , Dynamic
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Row | DenseMatrixOption::Vector> {
        using Base = Internal::AbstractDenseMatrixStorage<DenseMatrixStorage<T
                                                                           , DenseMatrixOption::Row | DenseMatrixOption::Vector
                                                                           , Dynamic
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                         , DenseMatrixOption::Row | DenseMatrixOption::Vector>;
    public:
        using Base::Base;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return Base::getLength(); }
        [[nodiscard]] size_t getColumn() const noexcept { return Base::operator[](0).getLength(); }
    };
}