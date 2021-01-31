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

#include "DenseMatrixStorageHelper.h"

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
            constexpr static int MatrixType = type;
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
    class DenseMatrixStorage<T, DenseMatrixType::Column | DenseMatrixType::Element, Row, Column, MaxRow, MaxColumn>
            : public Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Column | DenseMatrixType::Element
                                                                           , Row
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Column | DenseMatrixType::Element> {
        using Base = Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Column | DenseMatrixType::Element
                                                                           , Row
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Column | DenseMatrixType::Element>;
    protected:
        using Base::Base;
    public:
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return Row; }
        [[nodiscard]] constexpr static size_t getColumn() noexcept { return Column; }
    };

    template<class T, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixType::Column | DenseMatrixType::Element, Dynamic, Column, MaxRow, MaxColumn>
            : public Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Column | DenseMatrixType::Element
                                                                           , Dynamic
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Column | DenseMatrixType::Element> {
        using Base = Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Column | DenseMatrixType::Element
                                                                           , Dynamic
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Column | DenseMatrixType::Element>;
    protected:
        using Base::Base;
    public:
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { assert(Base::getLength() % getColumn() == 0); return Base::getLength() / getColumn(); }
        [[nodiscard]] constexpr static size_t getColumn() noexcept { return Column; }
    };

    template<class T, size_t Row, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixType::Column | DenseMatrixType::Element, Row, Dynamic, MaxRow, MaxColumn>
            : public Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Column | DenseMatrixType::Element
                                                                           , Row
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Column | DenseMatrixType::Element> {
        using Base = Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Column | DenseMatrixType::Element
                                                                           , Row
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Column | DenseMatrixType::Element>;
    protected:
        using Base::Base;
    public:
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return Row; }
        [[nodiscard]] size_t getColumn() const noexcept { assert(Base::getLength() % getRow() == 0); return Base::getLength() / getRow(); }
    };

    template<class T, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixType::Column | DenseMatrixType::Element, Dynamic, Dynamic, MaxRow, MaxColumn>
            : public Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Column | DenseMatrixType::Element
                                                                           , Dynamic
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Column | DenseMatrixType::Element> {
        using Base = Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Column | DenseMatrixType::Element
                                                                           , Dynamic
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Column | DenseMatrixType::Element>;
    protected:
        size_t row;
    public:
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return row; }
        [[nodiscard]] size_t getColumn() const noexcept { assert(Base::getLength() % getRow() == 0); return Base::getLength() / getRow(); }
    protected:
        DenseMatrixStorage() : Base(), row(0) {}
        DenseMatrixStorage(size_t row, size_t column, const T& t) : Base(row, column, t), row(0) {}
        DenseMatrixStorage(const DenseMatrixStorage&) = default;
        DenseMatrixStorage(DenseMatrixStorage&&) noexcept = default;
        ~DenseMatrixStorage() = default;
        /* Operators */
        DenseMatrixStorage& operator=(const DenseMatrixStorage&) = default;
        DenseMatrixStorage& operator=(DenseMatrixStorage&&) noexcept = default;
        /* Helper */
        void swap(DenseMatrixStorage& storage) noexcept { Base::swap(storage); }
    };
    ////////////////////////////////////RowElement////////////////////////////////////
    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixType::Row | DenseMatrixType::Element, Row, Column, MaxRow, MaxColumn>
            : public Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Row | DenseMatrixType::Element
                                                                           , Row
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Row | DenseMatrixType::Element> {
        using Base = Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Row | DenseMatrixType::Element
                                                                           , Row
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Row | DenseMatrixType::Element>;
    protected:
        using Base::Base;
    public:
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return Row; }
        [[nodiscard]] constexpr static size_t getColumn() noexcept { return Column; }
    };

    template<class T, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixType::Row | DenseMatrixType::Element, Dynamic, Column, MaxRow, MaxColumn>
            : public Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Row | DenseMatrixType::Element
                                                                           , Dynamic
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Row | DenseMatrixType::Element> {
        using Base = Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Row | DenseMatrixType::Element
                                                                           , Dynamic
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Row | DenseMatrixType::Element>;
    protected:
        using Base::Base;
    public:
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { assert(Base::getLength() % getColumn() == 0); return Base::getLength() / getColumn(); }
        [[nodiscard]] constexpr static size_t getColumn() noexcept { return Column; }
    };

    template<class T, size_t Row, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixType::Row | DenseMatrixType::Element, Row, Dynamic, MaxRow, MaxColumn>
            : public Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Row | DenseMatrixType::Element
                                                                           , Row
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Row | DenseMatrixType::Element> {
        using Base = Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Row | DenseMatrixType::Element
                                                                           , Row
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Row | DenseMatrixType::Element>;
    protected:
        using Base::Base;
    public:
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return Row; }
        [[nodiscard]] size_t getColumn() const noexcept { assert(Base::getLength() % getRow() == 0); return Base::getLength() / getRow(); }
    };

    template<class T, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixType::Row | DenseMatrixType::Element, Dynamic, Dynamic, MaxRow, MaxColumn>
            : public Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Row | DenseMatrixType::Element
                                                                           , Dynamic
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Row | DenseMatrixType::Element> {
        using Base = Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Row | DenseMatrixType::Element
                                                                           , Dynamic
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Row | DenseMatrixType::Element>;
    protected:
        size_t row;
    public:
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return row; }
        [[nodiscard]] size_t getColumn() const noexcept { assert(Base::getLength() % getRow() == 0); return Base::getLength() / getRow(); }
    protected:
        DenseMatrixStorage() : Base(), row(0) {}
        DenseMatrixStorage(size_t row, size_t column, const T& t) : Base(row, column, t), row(0) {}
        DenseMatrixStorage(const DenseMatrixStorage&) = default;
        DenseMatrixStorage(DenseMatrixStorage&&) noexcept = default;
        ~DenseMatrixStorage() = default;
        /* Operators */
        DenseMatrixStorage& operator=(const DenseMatrixStorage&) = default;
        DenseMatrixStorage& operator=(DenseMatrixStorage&&) noexcept = default;
        /* Helper */
        void swap(DenseMatrixStorage& storage) noexcept { Base::swap(storage); }
    };
    ////////////////////////////////////ColumnVector////////////////////////////////////
    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixType::Column | DenseMatrixType::Vector, Row, Column, MaxRow, MaxColumn>
            : public Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Column | DenseMatrixType::Vector
                                                                           , Row
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Column | DenseMatrixType::Vector> {
        using Base = Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Column | DenseMatrixType::Vector
                                                                           , Row
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Column | DenseMatrixType::Vector>;
    protected:
        using Base::Base;
    public:
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return Row; }
        [[nodiscard]] constexpr static size_t getColumn() noexcept { return Column; }
    };

    template<class T, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixType::Column | DenseMatrixType::Vector, Dynamic, Column, MaxRow, MaxColumn>
            : public Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Column | DenseMatrixType::Element
                                                                           , Dynamic
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Column | DenseMatrixType::Element> {
        using Base = Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Column | DenseMatrixType::Element
                                                                           , Dynamic
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Column | DenseMatrixType::Element>;
    protected:
        using Base::Base;
    public:
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return Base::operator[](0).getLength(); }
        [[nodiscard]] constexpr static size_t getColumn() noexcept { return Column; }
    };

    template<class T, size_t Row, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixType::Column | DenseMatrixType::Vector, Row, Dynamic, MaxRow, MaxColumn>
            : public Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Column | DenseMatrixType::Element
                                                                           , Row
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Column | DenseMatrixType::Element> {
        using Base = Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Column | DenseMatrixType::Element
                                                                           , Row
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Column | DenseMatrixType::Element>;
    protected:
        using Base::Base;
    public:
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return Row; }
        [[nodiscard]] size_t getColumn() const noexcept {return Base::getLength(); }
    };

    template<class T, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixType::Column | DenseMatrixType::Vector, Dynamic, Dynamic, MaxRow, MaxColumn>
            : public Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Column | DenseMatrixType::Vector
                                                                           , Dynamic
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Column | DenseMatrixType::Vector> {
        using Base = Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Column | DenseMatrixType::Vector
                                                                           , Dynamic
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Column | DenseMatrixType::Vector>;
    protected:
        using Base::Base;
    public:
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return Base::operator[](0).getLength(); }
        [[nodiscard]] size_t getColumn() const noexcept { return Base::getLength(); }
    };
    ////////////////////////////////////RowVector////////////////////////////////////
    template<class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixType::Row | DenseMatrixType::Vector, Row, Column, MaxRow, MaxColumn>
            : public Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Row | DenseMatrixType::Vector
                                                                           , Row
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Row | DenseMatrixType::Vector> {
        using Base = Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Row | DenseMatrixType::Vector
                                                                           , Row
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Row | DenseMatrixType::Vector>;
    protected:
        using Base::Base;
    public:
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return Row; }
        [[nodiscard]] constexpr static size_t getColumn() noexcept { return Column; }
    };

    template<class T, size_t Column, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixType::Row | DenseMatrixType::Vector, Dynamic, Column, MaxRow, MaxColumn>
            : public Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Row | DenseMatrixType::Element
                                                                           , Dynamic
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Row | DenseMatrixType::Element> {
        using Base = Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Row | DenseMatrixType::Element
                                                                           , Dynamic
                                                                           , Column
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Row | DenseMatrixType::Element>;
    protected:
        using Base::Base;
    public:
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return Base::getLength(); }
        [[nodiscard]] constexpr static size_t getColumn() noexcept { return Column; }
    };

    template<class T, size_t Row, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixType::Row | DenseMatrixType::Vector, Row, Dynamic, MaxRow, MaxColumn>
            : public Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Row | DenseMatrixType::Element
                                                                           , Row
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Row | DenseMatrixType::Element> {
        using Base = Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Row | DenseMatrixType::Element
                                                                           , Row
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Row | DenseMatrixType::Element>;
    protected:
        using Base::Base;
    public:
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return Row; }
        [[nodiscard]] size_t getColumn() const noexcept {return Base::operator[](0).getLength(); }
    };

    template<class T, size_t MaxRow, size_t MaxColumn>
    class DenseMatrixStorage<T, DenseMatrixType::Row | DenseMatrixType::Vector, Dynamic, Dynamic, MaxRow, MaxColumn>
            : public Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Row | DenseMatrixType::Vector
                                                                           , Dynamic
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Row | DenseMatrixType::Vector> {
        using Base = Internal::DenseMatrixStorageHelper<DenseMatrixStorage<T
                                                                           , DenseMatrixType::Row | DenseMatrixType::Vector
                                                                           , Dynamic
                                                                           , Dynamic
                                                                           , MaxRow
                                                                           , MaxColumn>
                                                        , DenseMatrixType::Row | DenseMatrixType::Vector>;
    protected:
        using Base::Base;
    public:
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return Base::getLength(); }
        [[nodiscard]] size_t getColumn() const noexcept { return Base::operator[](0).getLength(); }
    };
}