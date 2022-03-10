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

#include <cassert>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrixImpl/DenseMatrixOption.h"
#include "Physica/Utils/Template/CRTPBase.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

namespace Physica::Core::Internal {
    template<class T> class Traits;
    /**
     * This layer handles extends from \class Array, and decieds which function from \class Array
     * should expose to its child classes.
     */
    template<class T, size_t Size, size_t MaxSize>
    class DenseMatrixStorageHelper : private Utils::Array<T, Size, MaxSize> {
        static_assert(Size == MaxSize, "MaxSize that larger than Size makes no sence.");

        using Base = Utils::Array<T, Size, MaxSize>;
    public:
        using Base::Iterator;
        using Base::ConstIterator;
        using Base::ReverseIterator;
        using Base::ConstReverseIterator;
    public:
        using Base::Base;
        /* Iterator */
        using Base::begin;
        using Base::end;
        using Base::cbegin;
        using Base::cend;
        using Base::rbegin;
        using Base::rend;
        using Base::crbegin;
        using Base::crend;
        /* Operators */
        using Base::operator[];
    protected:
        /* Operations */
        using Base::resize;
        /* Getters */
        using Base::getLength;
        void swap(DenseMatrixStorageHelper& helper) noexcept { Base::swap(helper); }
    };

    template<class T, size_t MaxSize>
    class DenseMatrixStorageHelper<T, Dynamic, MaxSize> : private Utils::Array<T, Dynamic, MaxSize> {
        using Base = Utils::Array<T, Dynamic, MaxSize>;
    public:
        using Base::Iterator;
        using Base::ConstIterator;
        using Base::ReverseIterator;
        using Base::ConstReverseIterator;
    public:
        using Base::Base;
        /* Iterator */
        using Base::begin;
        using Base::end;
        using Base::cbegin;
        using Base::cend;
        using Base::rbegin;
        using Base::rend;
        using Base::crbegin;
        using Base::crend;
        /* Operators */
        using Base::operator[];
    protected:
        /* Operations */
        using Base::append;
        using Base::resize;
        /* Getters */
        using Base::getLength;
        void swap(DenseMatrixStorageHelper& helper) noexcept { Base::swap(helper); }
    };
    /**
     * This layer handles specialization of operator().
     */
    template<class Derived, int type>
    class AbstractDenseMatrixStorage;

    template<class Derived>
    class AbstractDenseMatrixStorage<Derived, DenseMatrixOption::Column | DenseMatrixOption::Element>
            : public DenseMatrixStorageHelper<typename Traits<Derived>::ScalarType
                                             , Traits<Derived>::SizeAtCompile
                                             , Traits<Derived>::MaxSizeAtCompile>
            , public Utils::CRTPBase<Derived> {
        static_assert(Traits<Derived>::MatrixOption == (DenseMatrixOption::Column | DenseMatrixOption::Element)
                      , "Invalid Derived type.");

        using Base = DenseMatrixStorageHelper<typename Traits<Derived>::ScalarType
                                             , Traits<Derived>::SizeAtCompile
                                             , Traits<Derived>::MaxSizeAtCompile>;
    private:
        using Utils::CRTPBase<Derived>::getDerived;
        using T = typename Traits<Derived>::ScalarType;
    public:
        AbstractDenseMatrixStorage() = default;
        AbstractDenseMatrixStorage(size_t row, size_t column) : Base(row * column) {}
        AbstractDenseMatrixStorage(size_t row, size_t column, const T& t) : Base(row * column, t) {}
        AbstractDenseMatrixStorage(std::initializer_list<T> list) : Base(list) {}
        /* Operators */
        [[nodiscard]] T& operator()(size_t r, size_t c) {
            assert(r < getDerived().getRow() && c < getDerived().getColumn());
            return Base::operator[](getDerived().getRow() * c + r);
        }

        [[nodiscard]] const T& operator()(size_t r, size_t c) const {
            assert(r < getDerived().getRow() && c < getDerived().getColumn());
            return Base::operator[](getDerived().getRow() * c + r);
        }
        /* Operations */
        void resize(size_t row, size_t column) { Base::resize(row * column); }
        template<size_t Length, size_t MaxLength>
        void appendRow(const Vector<T, Length, MaxLength>& v);
        void removeColumnAt(size_t index);
        void rowSwap(size_t r1, size_t r2);
        void columnSwap(size_t c1, size_t r1);
        /* Helpers */
        /**
         * Does not use Base::swap directly to avoid incorrect swaps
         */
        void swap(AbstractDenseMatrixStorage& storage) noexcept { Base::swap(storage); }
    };

    template<class Derived>
    class AbstractDenseMatrixStorage<Derived, DenseMatrixOption::Row | DenseMatrixOption::Element>
            : public DenseMatrixStorageHelper<typename Traits<Derived>::ScalarType
                                             , Traits<Derived>::SizeAtCompile
                                             , Traits<Derived>::MaxSizeAtCompile>
            , public Utils::CRTPBase<Derived> {
        static_assert(Traits<Derived>::MatrixOption == (DenseMatrixOption::Row | DenseMatrixOption::Element)
                      , "Invalid Derived type.");

        using Base = DenseMatrixStorageHelper<typename Traits<Derived>::ScalarType
                                             , Traits<Derived>::SizeAtCompile
                                             , Traits<Derived>::MaxSizeAtCompile>;
    private:
        using Utils::CRTPBase<Derived>::getDerived;
        using T = typename Traits<Derived>::ScalarType;
    public:
        AbstractDenseMatrixStorage() = default;
        AbstractDenseMatrixStorage(size_t row, size_t column) : Base(row * column) {}
        AbstractDenseMatrixStorage(size_t row, size_t column, const T& t) : Base(row * column, t) {}
        AbstractDenseMatrixStorage(std::initializer_list<T> list) : Base(list) {}
        /* Operators */
        [[nodiscard]] T& operator()(size_t r, size_t c) {
            assert(r < getDerived().getRow() && c < getDerived().getColumn());
            return Base::operator[](getDerived().getColumn() * r + c);
        }

        [[nodiscard]] const T& operator()(size_t r, size_t c) const  {
            assert(r < getDerived().getRow() && c < getDerived().getColumn());
            return Base::operator[](getDerived().getColumn() * r + c);
        }
        /* Operations */
        void resize(size_t row, size_t column) { Base::resize(row * column); }
        template<size_t Length, size_t MaxLength>
        void appendRow(const Vector<T, Length, MaxLength>& v);
        void removeColumnAt(size_t index);
        void rowSwap(size_t r1, size_t r2);
        void columnSwap(size_t c1, size_t r1);
        /* Helpers */
        void swap(AbstractDenseMatrixStorage& storage) noexcept { Base::swap(storage); }
    };

    template<class Derived>
    class AbstractDenseMatrixStorage<Derived, DenseMatrixOption::Column | DenseMatrixOption::Vector>
            : public DenseMatrixStorageHelper<Vector<typename Traits<Derived>::ScalarType, Traits<Derived>::RowAtCompile, Traits<Derived>::MaxRowAtCompile>
                                             , Traits<Derived>::ColumnAtCompile
                                             , Traits<Derived>::MaxColumnAtCompile>
            , public Utils::CRTPBase<Derived> {
        static_assert(Traits<Derived>::MatrixOption == (DenseMatrixOption::Column | DenseMatrixOption::Vector)
                      , "Invalid Derived type.");
    public:
        using VectorType = Vector<typename Traits<Derived>::ScalarType, Traits<Derived>::RowAtCompile, Traits<Derived>::MaxRowAtCompile>;
    private:
        using Base = DenseMatrixStorageHelper<VectorType
                                             , Traits<Derived>::ColumnAtCompile
                                             , Traits<Derived>::MaxColumnAtCompile>;
    private:
        using Utils::CRTPBase<Derived>::getDerived;
        using T = typename Traits<Derived>::ScalarType;
    public:
        AbstractDenseMatrixStorage() = default;
        AbstractDenseMatrixStorage(size_t row, size_t column) : Base(column, row) {}
        AbstractDenseMatrixStorage(size_t row, size_t column, const T& t) : Base(column, row, t) {}
        AbstractDenseMatrixStorage(std::initializer_list<VectorType> list) : Base(list) {}
        /* Operators */
        [[nodiscard]] T& operator()(size_t r, size_t c) {
            assert(r < getDerived().getRow() && c < getDerived().getColumn());
            return Base::operator[](c)[r];
        }

        [[nodiscard]] const T& operator()(size_t r, size_t c) const  {
            assert(r < getDerived().getRow() && c < getDerived().getColumn());
            return Base::operator[](c)[r];
        }
        /* Operations */
        void resize(size_t row, size_t column);
        template<size_t Length, size_t MaxLength>
        void appendRow(const Vector<T, Length, MaxLength>& v);
        void removeColumnAt(size_t index);
        void rowSwap(size_t r1, size_t r2);
        void columnSwap(size_t c1, size_t r1);
        /* Helpers */
        void swap(AbstractDenseMatrixStorage& storage) noexcept { Base::swap(storage); }
    };

    template<class Derived>
    class AbstractDenseMatrixStorage<Derived, DenseMatrixOption::Row | DenseMatrixOption::Vector>
            : public DenseMatrixStorageHelper<Vector<typename Traits<Derived>::ScalarType, Traits<Derived>::ColumnAtCompile, Traits<Derived>::MaxColumnAtCompile>
                                             , Traits<Derived>::RowAtCompile
                                             , Traits<Derived>::MaxRowAtCompile>
            , public Utils::CRTPBase<Derived> {
        static_assert(Traits<Derived>::MatrixOption == (DenseMatrixOption::Row | DenseMatrixOption::Vector)
                      , "Invalid Derived type.");
    public:
        using VectorType = Vector<typename Traits<Derived>::ScalarType, Traits<Derived>::ColumnAtCompile, Traits<Derived>::MaxColumnAtCompile>;
    private:
        using Base = DenseMatrixStorageHelper<VectorType
                                             , Traits<Derived>::RowAtCompile
                                             , Traits<Derived>::MaxRowAtCompile>;
    private:
        using Utils::CRTPBase<Derived>::getDerived;
        using T = typename Traits<Derived>::ScalarType;
    public:
        AbstractDenseMatrixStorage() = default;
        AbstractDenseMatrixStorage(size_t row, size_t column) : Base(row, column) {}
        AbstractDenseMatrixStorage(size_t row, size_t column, const T& t) : Base(row, column, t) {}
        AbstractDenseMatrixStorage(std::initializer_list<VectorType> list) : Base(list) {}
        /* Operators */
        [[nodiscard]] T& operator()(size_t r, size_t c) {
            assert(r < getDerived().getRow() && c < getDerived().getColumn());
            return Base::operator[](r)[c];
        }

        [[nodiscard]] const T& operator()(size_t r, size_t c) const  {
            assert(r < getDerived().getRow() && c < getDerived().getColumn());
            return Base::operator[](r)[c];
        }
        /* Operations */
        void resize(size_t row, size_t column);
        template<size_t Length, size_t MaxLength>
        void appendRow(const Vector<T, Length, MaxLength>& v);
        void removeColumnAt(size_t index);
        void rowSwap(size_t r1, size_t r2);
        void columnSwap(size_t c1, size_t r1);
        /* Helpers */
        void swap(AbstractDenseMatrixStorage& storage) noexcept { Base::swap(storage); }
    };
}

#include "AbstractDenseMatrixStorageImpl.h"
