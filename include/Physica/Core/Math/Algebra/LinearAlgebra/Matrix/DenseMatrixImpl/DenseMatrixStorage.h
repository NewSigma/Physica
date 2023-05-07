/*
 * Copyright 2021-2022 WeiBo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixOption.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

namespace Physica::Core::Internal {
    template<class T> class Traits;
    /**
     * This layer handles specialization of operator().
     */
    template<class Derived, int type>
    class DenseMatrixStorage;

    template<class Derived>
    class DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Element>
            : public Utils::Array<typename Traits<Derived>::ScalarType,
                                  Traits<Derived>::SizeAtCompile,
                                  Traits<Derived>::MaxSizeAtCompile>
            , public Utils::CRTPBase<Derived, 1> {
        static_assert(Traits<Derived>::Option == (MatrixOption::Column | MatrixOption::Element)
                      , "Invalid Derived type.");
        using This = DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Element>;
    public:
        using Base = Utils::Array<typename Traits<Derived>::ScalarType, Traits<Derived>::SizeAtCompile, Traits<Derived>::MaxSizeAtCompile>;
    private:
        using Utils::CRTPBase<Derived, 1>::getDerived;
        using T = typename Traits<Derived>::ScalarType;
    public:
        DenseMatrixStorage() = default;
        DenseMatrixStorage(size_t row, size_t column) : Base(row * column) {}
        DenseMatrixStorage(size_t row, size_t column, const T& t) : Base(row * column, t) {}
        DenseMatrixStorage(std::initializer_list<T> list) : Base(list) {}
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
        [[nodiscard]] device_obj<This> toDevice() const;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return Base::getLength(); }
    private:
        DenseMatrixStorage(Base array) : Base(std::move(array)) {}

        using Base::getCapacity;
        friend class device_obj<This>;
    };

    template<class Derived>
    class DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Element>
            : public Utils::Array<typename Traits<Derived>::ScalarType,
                                  Traits<Derived>::SizeAtCompile,
                                  Traits<Derived>::MaxSizeAtCompile>
            , public Utils::CRTPBase<Derived, 1> {
        static_assert(Traits<Derived>::Option == (MatrixOption::Row | MatrixOption::Element)
                      , "Invalid Derived type.");
        using This = DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Element>;
    public:
        using Base = Utils::Array<typename Traits<Derived>::ScalarType, Traits<Derived>::SizeAtCompile, Traits<Derived>::MaxSizeAtCompile>;
    private:
        using Utils::CRTPBase<Derived, 1>::getDerived;
        using T = typename Traits<Derived>::ScalarType;
    public:
        DenseMatrixStorage() = default;
        DenseMatrixStorage(size_t row, size_t column) : Base(row * column) {}
        DenseMatrixStorage(size_t row, size_t column, const T& t) : Base(row * column, t) {}
        DenseMatrixStorage(std::initializer_list<T> list) : Base(list) {}
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
        [[nodiscard]] device_obj<This> toDevice() const;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return Base::getLength(); }
    private:
        DenseMatrixStorage(Base array) : Base(std::move(array)) {}

        using Base::getCapacity;
        friend class device_obj<This>;
    };

    template<class Derived>
    class DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Vector>
            : public Utils::Array<Vector<typename Traits<Derived>::ScalarType, Traits<Derived>::RowAtCompile, Traits<Derived>::MaxRowAtCompile>,
                                  Traits<Derived>::ColumnAtCompile,
                                  Traits<Derived>::MaxColumnAtCompile>
            , public Utils::CRTPBase<Derived, 1> {
        static_assert(Traits<Derived>::Option == (MatrixOption::Column | MatrixOption::Vector)
                      , "Invalid Derived type.");
        using This = DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Vector>;
    public:
        using VectorType = Vector<typename Traits<Derived>::ScalarType, Traits<Derived>::RowAtCompile, Traits<Derived>::MaxRowAtCompile>;
        using Base = Utils::Array<VectorType, Traits<Derived>::ColumnAtCompile, Traits<Derived>::MaxColumnAtCompile>;
        using Utils::CRTPBase<Derived, 1>::getDerived;
    private:
        using T = typename Traits<Derived>::ScalarType;
    public:
        DenseMatrixStorage() = default;
        DenseMatrixStorage(size_t row, size_t column) : Base(column, row) {}
        DenseMatrixStorage(size_t row, size_t column, const T& t) : Base(column, row, t) {}
        DenseMatrixStorage(std::initializer_list<VectorType> list) : Base(list) {}
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
        [[nodiscard]] device_obj<This> toDevice() const;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept {
            return Base::getLength() == 0 ? 0 : Base::getLength() * Base::operator[](0).getLength();
        }
    private:
        DenseMatrixStorage(Base array) : Base(std::move(array)) {}

        using Base::getCapacity;
        friend class device_obj<This>;
    };

    template<class Derived>
    class DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>
            : public Utils::Array<Vector<typename Traits<Derived>::ScalarType, Traits<Derived>::ColumnAtCompile, Traits<Derived>::MaxColumnAtCompile>,
                                         Traits<Derived>::RowAtCompile,
                                         Traits<Derived>::MaxRowAtCompile>
            , public Utils::CRTPBase<Derived, 1> {
        static_assert(Traits<Derived>::Option == (MatrixOption::Row | MatrixOption::Vector)
                      , "Invalid Derived type.");
        using This = DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>;
    public:
        using VectorType = Vector<typename Traits<Derived>::ScalarType, Traits<Derived>::ColumnAtCompile, Traits<Derived>::MaxColumnAtCompile>;
        using Base = Utils::Array<VectorType, Traits<Derived>::RowAtCompile, Traits<Derived>::MaxRowAtCompile>;
    private:
        using Utils::CRTPBase<Derived, 1>::getDerived;
        using T = typename Traits<Derived>::ScalarType;
    public:
        DenseMatrixStorage() = default;
        DenseMatrixStorage(size_t row, size_t column) : Base(row, column) {}
        DenseMatrixStorage(size_t row, size_t column, const T& t) : Base(row, column, t) {}
        DenseMatrixStorage(std::initializer_list<VectorType> list) : Base(list) {}
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
        [[nodiscard]] device_obj<This> toDevice() const;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept {
            return Base::getLength() == 0 ? 0 : Base::getLength() * Base::operator[](0).getLength();
        }
    private:
        DenseMatrixStorage(Base array) : Base(std::move(array)) {}

        using Base::getCapacity;
        friend class device_obj<This>;
    };
}

#include "DenseMatrixStorageImpl.h"
