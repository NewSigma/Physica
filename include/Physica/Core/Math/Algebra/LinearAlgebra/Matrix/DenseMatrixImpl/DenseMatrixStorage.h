/*
 * Copyright 2021-2023 WeiBo He.
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
                                  Traits<Derived>::MaxSizeAtCompile,
                                  typename Traits<Derived>::AllocatorType>
            , public Utils::CRTPBase<Derived, 1> {
        static_assert(Traits<Derived>::Option == (MatrixOption::Column | MatrixOption::Element)
                      , "Invalid Derived type.");
        using This = DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Element>;
        using T = typename Traits<Derived>::ScalarType;
    public:
        using Base = Utils::Array<typename Traits<Derived>::ScalarType, Traits<Derived>::SizeAtCompile, Traits<Derived>::MaxSizeAtCompile, typename Traits<Derived>::AllocatorType>;
        using InitializerType = T;
    private:
        using Utils::CRTPBase<Derived, 1>::getDerived;
    public:
        DenseMatrixStorage() = default;
        DenseMatrixStorage(size_t row, size_t column) : Base(row * column) {}
        DenseMatrixStorage(size_t row, size_t column, T value) : Base(row * column, value) {}
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
                                  Traits<Derived>::MaxSizeAtCompile,
                                  typename Traits<Derived>::AllocatorType>
            , public Utils::CRTPBase<Derived, 1> {
        static_assert(Traits<Derived>::Option == (MatrixOption::Row | MatrixOption::Element)
                      , "Invalid Derived type.");
        using This = DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Element>;
        using T = typename Traits<Derived>::ScalarType;
    public:
        using Base = Utils::Array<typename Traits<Derived>::ScalarType, Traits<Derived>::SizeAtCompile, Traits<Derived>::MaxSizeAtCompile, typename Traits<Derived>::AllocatorType>;
        using InitializerType = T;
    private:
        using Utils::CRTPBase<Derived, 1>::getDerived;
    public:
        DenseMatrixStorage() = default;
        DenseMatrixStorage(size_t row, size_t column) : Base(row * column) {}
        DenseMatrixStorage(size_t row, size_t column, T value) : Base(row * column, value) {}
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
            : public Utils::CRTPBase<Derived, 1> {
        static_assert(Traits<Derived>::Option == (MatrixOption::Column | MatrixOption::Vector)
                      , "Invalid Derived type.");
        using This = DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Vector>;
        using Base = Utils::CRTPBase<Derived, 1>;
        using T = typename Traits<Derived>::ScalarType;
        using AllocatorTypeT = typename Traits<Derived>::AllocatorType;
    public:
        using VectorType = Vector<typename Traits<Derived>::ScalarType, Traits<Derived>::RowAtCompile, Traits<Derived>::MaxRowAtCompile, AllocatorTypeT>;
        using InitializerType = VectorType;
    private:
        using AllocatorTypeV = typename Utils::ChangeAllocatorValueType<AllocatorTypeT, VectorType>::Type;
        using ArrayType = Utils::Array<VectorType, Traits<Derived>::ColumnAtCompile, Traits<Derived>::MaxColumnAtCompile, AllocatorTypeV>;

        ArrayType array;
    public:
        DenseMatrixStorage() = default;
        DenseMatrixStorage(size_t row, size_t column) : array(column, row) {}
        DenseMatrixStorage(size_t row, size_t column, T value) : array(column, row, value) {}
        DenseMatrixStorage(std::initializer_list<VectorType> list) : array(list) {}
        /* Operators */
        [[nodiscard]] T& operator()(size_t r, size_t c) {
            assert(r < Base::getDerived().getRow() && c < Base::getDerived().getColumn());
            return array[c][r];
        }

        [[nodiscard]] const T& operator()(size_t r, size_t c) const  {
            assert(r < Base::getDerived().getRow() && c < Base::getDerived().getColumn());
            return array[c][r];
        }
        /* Operations */
        void resize(size_t row, size_t column);
        template<size_t Length, size_t MaxLength>
        void appendRow(const Vector<T, Length, MaxLength>& v);
        void rowSwap(size_t r1, size_t r2);
        void columnSwap(size_t c1, size_t r1);
        [[nodiscard]] device_obj<This> toDevice() const;
        void swap(DenseMatrixStorage& obj) noexcept { array.swap(obj.array); }
        /* Getters */
        [[nodiscard]] ArrayType& asArray() noexcept { return array; }
        [[nodiscard]] const ArrayType& asArray() const noexcept { return array; }
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept {
            return array.getLength() == 0 ? 0 : array.getLength() * array[0].getLength();
        }
    private:
        DenseMatrixStorage(ArrayType array_) : array(std::move(array_)) {}

        friend class device_obj<This>;
    };

    template<class Derived>
    class DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>
            : public Utils::CRTPBase<Derived, 1> {
        static_assert(Traits<Derived>::Option == (MatrixOption::Row | MatrixOption::Vector)
                      , "Invalid Derived type.");
        using This = DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>;
        using Base = Utils::CRTPBase<Derived, 1>;
        using T = typename Traits<Derived>::ScalarType;
        using AllocatorTypeT = typename Traits<Derived>::AllocatorType;
    public:
        using VectorType = Vector<typename Traits<Derived>::ScalarType, Traits<Derived>::ColumnAtCompile, Traits<Derived>::MaxColumnAtCompile, AllocatorTypeT>;
        using InitializerType = VectorType;
    private:
        using AllocatorTypeV = typename Utils::ChangeAllocatorValueType<AllocatorTypeT, VectorType>::Type;
        using ArrayType = Utils::Array<VectorType, Traits<Derived>::RowAtCompile, Traits<Derived>::MaxRowAtCompile, AllocatorTypeV>;
        
        ArrayType array;
    public:
        DenseMatrixStorage() = default;
        DenseMatrixStorage(size_t row, size_t column) : array(row, column) {}
        DenseMatrixStorage(size_t row, size_t column, T value) : array(row, column, value) {}
        DenseMatrixStorage(std::initializer_list<VectorType> list) : array(list) {}
        /* Operators */
        [[nodiscard]] T& operator()(size_t r, size_t c) {
            assert(r < Base::getDerived().getRow() && c < Base::getDerived().getColumn());
            return array[r][c];
        }

        [[nodiscard]] const T& operator()(size_t r, size_t c) const  {
            assert(r < Base::getDerived().getRow() && c < Base::getDerived().getColumn());
            return array[r][c];
        }
        /* Operations */
        void resize(size_t row, size_t column);
        template<size_t Length, size_t MaxLength>
        void appendRow(const Vector<T, Length, MaxLength>& v);
        void rowSwap(size_t r1, size_t r2);
        void columnSwap(size_t c1, size_t r1);
        [[nodiscard]] device_obj<This> toDevice() const;
        void swap(DenseMatrixStorage& obj) noexcept { array.swap(obj.array); }
        /* Getters */
        [[nodiscard]] ArrayType& asArray() noexcept { return array; }
        [[nodiscard]] const ArrayType& asArray() const noexcept { return array; }
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept {
            return array.getLength() == 0 ? 0 : array.getLength() * array[0].getLength();
        }
    private:
        DenseMatrixStorage(ArrayType array_) : array(std::move(array_)) {}

        friend class device_obj<This>;
    };
}

#include "DenseMatrixStorageImpl.h"
