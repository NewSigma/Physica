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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrixImpl/DenseMatrixType.h"
#include "Physica/Utils/Template/CRTPBase.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector.h"

namespace Physica::Core::Internal {
    template<class T> class Traits;
    /**
     * This layer handles extends from \class Array, and decieds which function from \class Array
     * should expose to its child classes.
     */
    template<class T, size_t Size, size_t MaxSize>
    class AbstractDenseMatrixStorage : private Utils::Array<T, Size, MaxSize> {
        using Base = Utils::Array<T, Size, MaxSize>;
    protected:
        using Base::Base;
    public:
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
        AbstractDenseMatrixStorage() = default;
        AbstractDenseMatrixStorage(const AbstractDenseMatrixStorage&) = default;
        AbstractDenseMatrixStorage(AbstractDenseMatrixStorage&&) = default;
        ~AbstractDenseMatrixStorage() = default;
        /* Operators */
        AbstractDenseMatrixStorage& operator=(const AbstractDenseMatrixStorage&) = default;
        AbstractDenseMatrixStorage& operator=(AbstractDenseMatrixStorage&& m) noexcept = default;
        /* Operations */
        using Base::append;
        /* Getters */
        using Base::getLength;
    };
    /**
     * This layre handles specialization of operator().
     */
    template<class Derived, int type>
    class DenseMatrixStorageHelper;

    template<class Derived>
    class DenseMatrixStorageHelper<Derived, DenseMatrixType::Column | DenseMatrixType::Element>
            : public AbstractDenseMatrixStorage<typename Traits<Derived>::ScalarType
                                                , Traits<Derived>::SizeAtCompile
                                                , Traits<Derived>::MaxSizeAtCompile>
            , public Utils::CRTPBase<Derived> {
        static_assert(Traits<Derived>::MatrixType == (DenseMatrixType::Column | DenseMatrixType::Element)
                      , "Invalid Derived type.");
    private:
        using Base = AbstractDenseMatrixStorage<typename Traits<Derived>::ScalarType
                                                , Traits<Derived>::SizeAtCompile
                                                , Traits<Derived>::MaxSizeAtCompile>;
        using Utils::CRTPBase<Derived>::getDerived;
        using T = typename Traits<Derived>::ScalarType;
    protected:
        using Base::Base;
    public:
        /* Operators */
        [[nodiscard]] T& operator()(size_t r, size_t c) {
            assert(r < getDerived().getRow() && c < getDerived().getColumn());
            return Base::operator[](getDerived().getRow() * r + c);
        }

        [[nodiscard]] const T& operator()(size_t r, size_t c) const  {
            assert(r < getDerived().getRow() && c < getDerived().getColumn());
            return Base::operator[](getDerived().getRow() * r + c);
        }
        /* Operations */
        template<size_t Length, size_t MaxLength>
        void appendRow(const Vector<T, Length, MaxLength>& v);
        void removeColumnAt(size_t index);
        void rowSwap(size_t r1, size_t r2);
    protected:
        DenseMatrixStorageHelper(size_t row, size_t column, const T& t) : Base(row * column, t) {}
        /* Helpers */
        void swap(DenseMatrixStorageHelper& helper) noexcept { Base::swap(helper); }
    };

    template<class Derived>
    class DenseMatrixStorageHelper<Derived, DenseMatrixType::Row | DenseMatrixType::Element>
            : public AbstractDenseMatrixStorage<typename Traits<Derived>::ScalarType
                                                , Traits<Derived>::SizeAtCompile
                                                , Traits<Derived>::MaxSizeAtCompile>
            , public Utils::CRTPBase<Derived> {
        static_assert(Traits<Derived>::MatrixType == (DenseMatrixType::Row | DenseMatrixType::Element)
                      , "Invalid Derived type.");
    private:
        using Base = AbstractDenseMatrixStorage<typename Traits<Derived>::ScalarType
                                                , Traits<Derived>::SizeAtCompile
                                                , Traits<Derived>::MaxSizeAtCompile>;
        using Utils::CRTPBase<Derived>::getDerived;
        using T = typename Traits<Derived>::ScalarType;
    protected:
        using Base::Base;
    public:
        /* Operators */
        [[nodiscard]] T& operator()(size_t r, size_t c) {
            assert(r < getDerived().getRow() && c < getDerived().getColumn());
            return Base::operator[](getDerived().getColumn() * c + r);
        }

        [[nodiscard]] const T& operator()(size_t r, size_t c) const  {
            assert(r < getDerived().getRow() && c < getDerived().getColumn());
            return Base::operator[](getDerived().getColumn() * c + r);
        }
        /* Operations */
        template<size_t Length, size_t MaxLength>
        void appendRow(const Vector<T, Length, MaxLength>& v);
        void removeColumnAt(size_t index);
        void rowSwap(size_t r1, size_t r2);
    protected:
        DenseMatrixStorageHelper(size_t row, size_t column, const T& t) : Base(row * column, t) {}
        /* Helpers */
        void swap(DenseMatrixStorageHelper& helper) noexcept { Base::swap(helper); }
    };

    template<class Derived>
    class DenseMatrixStorageHelper<Derived, DenseMatrixType::Column | DenseMatrixType::Vector>
            : public AbstractDenseMatrixStorage<Vector<typename Traits<Derived>::ScalarType, Traits<Derived>::RowAtCompile, Traits<Derived>::MaxRowAtCompile>
                                                , Traits<Derived>::ColumnAtCompile
                                                , Traits<Derived>::MaxColumnAtCompile>
            , public Utils::CRTPBase<Derived> {
        static_assert(Traits<Derived>::MatrixType == (DenseMatrixType::Column | DenseMatrixType::Vector)
                      , "Invalid Derived type.");
    public:
        using VectorType = Vector<typename Traits<Derived>::ScalarType, Traits<Derived>::RowAtCompile, Traits<Derived>::MaxRowAtCompile>;
    private:
        using Base = AbstractDenseMatrixStorage<VectorType
                                                , Traits<Derived>::ColumnAtCompile
                                                , Traits<Derived>::MaxColumnAtCompile>;
        using Utils::CRTPBase<Derived>::getDerived;
        using T = typename Traits<Derived>::ScalarType;
    protected:
        using Base::Base;
    public:
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
        template<size_t Length, size_t MaxLength>
        void appendRow(const Vector<T, Length, MaxLength>& v);
        void removeColumnAt(size_t index);
        void rowSwap(size_t r1, size_t r2);
    protected:
        DenseMatrixStorageHelper(size_t row, size_t column, const T& t) : Base(column, VectorType(row, t)) {}
        DenseMatrixStorageHelper(std::initializer_list<VectorType> list) : Base(std::move(list)) {}
        /* Helpers */
        void swap(DenseMatrixStorageHelper& helper) noexcept { Base::swap(helper); }
    };

    template<class Derived>
    class DenseMatrixStorageHelper<Derived, DenseMatrixType::Row | DenseMatrixType::Vector>
            : public AbstractDenseMatrixStorage<Vector<typename Traits<Derived>::ScalarType, Traits<Derived>::ColumnAtCompile, Traits<Derived>::MaxColumnAtCompile>
                                                , Traits<Derived>::RowAtCompile
                                                , Traits<Derived>::MaxRowAtCompile>
            , public Utils::CRTPBase<Derived> {
        static_assert(Traits<Derived>::MatrixType == (DenseMatrixType::Row | DenseMatrixType::Vector)
                      , "Invalid Derived type.");
    public:
        using VectorType = Vector<typename Traits<Derived>::ScalarType, Traits<Derived>::ColumnAtCompile, Traits<Derived>::MaxColumnAtCompile>;
    private:
        using Base = AbstractDenseMatrixStorage<VectorType
                                                , Traits<Derived>::RowAtCompile
                                                , Traits<Derived>::MaxRowAtCompile>;
        using Utils::CRTPBase<Derived>::getDerived;
        using T = typename Traits<Derived>::ScalarType;
    protected:
        using Base::Base;
    public:
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
        template<size_t Length, size_t MaxLength>
        void appendRow(const Vector<T, Length, MaxLength>& v);
        void removeColumnAt(size_t index);
        void rowSwap(size_t r1, size_t r2);
    protected:
        DenseMatrixStorageHelper(size_t row, size_t column, const T& t) : Base(row, VectorType(column, t)) {}
        DenseMatrixStorageHelper(std::initializer_list<VectorType> list) : Base(std::move(list)) {}
        /* Helpers */
        void swap(DenseMatrixStorageHelper& helper) noexcept { Base::swap(helper); }
    };
}

#include "DenseMatrixStorageHelperImpl.h"
