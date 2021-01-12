/*
 * Copyright 2020-2021 WeiBo He.
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
#ifndef PHYSICA_MATRIX_H
#define PHYSICA_MATRIX_H

#include <memory>
#include "MatrixBase.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector.h"

namespace Physica::Core {
    using Utils::Dynamic;
    /**
     * This enum decides how the data is stored in a matrix.
     * A dense matrix can be stored as elements or vectors in rows or columns.
     *
     * It is recommended that store elements to store a small matrix.
     */
    class DenseMatrixType {
    public:
        enum {
            Column = 0b00,
            Row = 0b01,
            Vector = 0b00,
            Element = 0b10
        };
    };
    /**
     * DenseMatrix class
     * A matrix can be either fixed matrix, which have its max size defined,
     * or dynamic matrix, whose size is dynamically changed.
     * 
     * \tparam type
     * type is combinations of \enum DenseMatrixType
     */
    template<class T = MultiScalar, int type = DenseMatrixType::Column | DenseMatrixType::Vector
            , size_t maxRow = Dynamic, size_t maxColumn = Dynamic>
    class DenseMatrix {
        static_assert(type == 0 || type == 1 || type == 2 || type == 3, "Invalid dense matrix type.");
    };
    /*!
     * Specialization for ElementColumn matrix.
     */
    template<class T, size_t maxRow, size_t maxColumn>
    class DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Element, maxRow, maxColumn>
            : private Utils::CStyleArray<T, maxRow * maxColumn> {
    public:
        using Base = Utils::CStyleArray<T, maxRow * maxColumn>;
    private:
        size_t row;
        size_t column;
    public:
        DenseMatrix() : row(0), column(0) {}
        DenseMatrix(size_t row, size_t column);
        DenseMatrix(std::initializer_list<Utils::CStyleArray<T, maxRow>> list);
        DenseMatrix(const DenseMatrix& matrix) = default;
        DenseMatrix(DenseMatrix&& matrix) noexcept;
        ~DenseMatrix() = default;
        /* Operators */
        [[nodiscard]] T& operator()(size_t r, size_t c) { return Base::operator[](row * r + c); }
        [[nodiscard]] const T& operator()(size_t r, size_t c) const  { return Base::operator[](row * r + c); }
        DenseMatrix& operator=(const DenseMatrix& m) = default;
        DenseMatrix& operator=(DenseMatrix&& m) noexcept;
        /* Iterators */
        typename Base::Iterator begin() { return Base::begin(); } //NOLINT virtual is not essential.
        typename Base::Iterator end() { return Base::end(); } //NOLINT virtual is not essential.
        /* Getters */
        [[nodiscard]] constexpr static int getType() { return DenseMatrixType::Column | DenseMatrixType::Element; }
        [[nodiscard]] size_t getRow() const { return row; }
        [[nodiscard]] size_t getColumn() const { return column; }
    };
    /*!
     * Specialization for ElementRow matrix.
     */
    template<class T, size_t maxRow, size_t maxColumn>
    class DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Element, maxRow, maxColumn>
            : private Utils::CStyleArray<T, maxRow * maxColumn> {
    public:
        using Base = Utils::CStyleArray<T, maxRow * maxColumn>;
    private:
        size_t row;
        size_t column;
    public:
        DenseMatrix() : row(0), column(0) {}
        DenseMatrix(size_t row, size_t column);
        DenseMatrix(std::initializer_list<Utils::CStyleArray<T, maxColumn>> list);
        DenseMatrix(const DenseMatrix& matrix) = default;
        DenseMatrix(DenseMatrix&& matrix) noexcept;
        ~DenseMatrix() = default;
        /* Operators */
        [[nodiscard]] T& operator()(size_t r, size_t c) { return Base::operator[](column * c + r); }
        [[nodiscard]] const T& operator()(size_t r, size_t c) const  { return Base::operator[](column * c + r); }
        DenseMatrix& operator=(const DenseMatrix& m) = default;
        DenseMatrix& operator=(DenseMatrix&& m) noexcept;
        /* Iterators */
        typename Base::Iterator begin() { return Base::begin(); } //NOLINT virtual is not essential.
        typename Base::Iterator end() { return Base::end(); } //NOLINT virtual is not essential.
        /* Getters */
        [[nodiscard]] constexpr static int getType() { return DenseMatrixType::Row | DenseMatrixType::Element; }
        [[nodiscard]] size_t getRow() const { return row; }
        [[nodiscard]] size_t getColumn() const { return column; }
    };
    /*!
     * Specialization for VectorColumn matrix.
     */
    template<class T, size_t maxRow, size_t maxColumn>
    class DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, maxColumn>
            : public Utils::CStyleArray<Vector<T, maxRow>, maxColumn> {
    public:
        using VectorType = Vector<T, maxRow>;
        using Base = Utils::CStyleArray<VectorType, maxColumn>;
    public:
        DenseMatrix() = default;
        explicit DenseMatrix(size_t length);
        DenseMatrix(std::initializer_list<VectorType> list);
        DenseMatrix(const DenseMatrix& matrix) = default;
        DenseMatrix(DenseMatrix&& matrix) noexcept;
        ~DenseMatrix() = default;
        /* Operators */
        [[nodiscard]] T& operator()(size_t r, size_t c) { return Base::operator[](c)[r]; }
        [[nodiscard]] const T& operator()(size_t r, size_t c) const  { return Base::operator[](c)[r]; }
        DenseMatrix& operator=(const DenseMatrix& m) = default;
        DenseMatrix& operator=(DenseMatrix&& m) noexcept;
        /* Iterators */
        typename Base::Iterator begin() { return Base::begin(); } //NOLINT virtual is not essential.
        typename Base::Iterator end() { return Base::end(); } //NOLINT virtual is not essential.
        /* Matrix Operations */
        void appendRow(const VectorType& v);
        void appendRow(VectorType&& v) noexcept;
        void appendColumn(const VectorType& v);
        void appendColumn(VectorType&& v) noexcept;
        void appendMatrixRow(const DenseMatrix& m);
        void appendMatrixRow(DenseMatrix&& m);
        void appendMatrixColumn(const DenseMatrix& m);
        void appendMatrixColumn(DenseMatrix&& m);
        void removeRowAt(size_t index);
        inline void removeColumnAt(size_t index);
        VectorType cutRow();
        VectorType cutColumn();
        DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, Dynamic, maxColumn> cutMatrixRow(size_t from);
        DenseMatrix<T, DenseMatrixType::Column | DenseMatrixType::Vector, maxRow, Dynamic> cutMatrixColumn(size_t from);
        void rowSwap(size_t r1, size_t r2) noexcept;
        void columnSwap(size_t c1, size_t c2) noexcept;
        void rowReduce(size_t r1, size_t r2, size_t element);
        void columnReduce(size_t c1, size_t c2, size_t element);
        /* Getters */
        [[nodiscard]] constexpr static int getType() { return DenseMatrixType::Column | DenseMatrixType::Vector; }
        [[nodiscard]] size_t getRow() const { return Base::operator[](0).getLength(); }
        [[nodiscard]] size_t getColumn() const { return Base::getLength(); }
    };
    /*!
     * Specialization for VectorRow matrix.
     */
    template<class T, size_t maxRow, size_t maxColumn>
    class DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, maxColumn>
            : public Utils::CStyleArray<Vector<T, maxColumn>, maxRow> {
    public:
        using VectorType = Vector<T, maxColumn>;
        using Base = Utils::CStyleArray<VectorType, maxRow>;
    public:
        DenseMatrix() = default;
        explicit DenseMatrix(size_t length);
        DenseMatrix(std::initializer_list<VectorType> list);
        DenseMatrix(const DenseMatrix& matrix) = default;
        DenseMatrix(DenseMatrix&& matrix) noexcept;
        ~DenseMatrix() = default;
        /* Operators */
        [[nodiscard]] T& operator()(size_t r, size_t c) { return Base::operator[](r)[c]; }
        [[nodiscard]] const T& operator()(size_t r, size_t c) const  { return Base::operator[](r)[c]; }
        DenseMatrix& operator=(const DenseMatrix& m) = default;
        DenseMatrix& operator=(DenseMatrix&& m) noexcept;
        /* Iterators */
        typename Base::Iterator begin() { return Base::begin(); } //NOLINT virtual is not essential.
        typename Base::Iterator end() { return Base::end(); } //NOLINT virtual is not essential.
        /* Matrix Operations */
        void appendRow(const VectorType& v);
        void appendRow(VectorType&& v) noexcept;
        void appendColumn(const VectorType& v);
        void appendColumn(VectorType&& v) noexcept;
        void appendMatrixRow(const DenseMatrix& m);
        void appendMatrixRow(DenseMatrix&& m);
        void appendMatrixColumn(const DenseMatrix& m);
        void appendMatrixColumn(DenseMatrix&& m);
        inline void removeRowAt(size_t index);
        void removeColumnAt(size_t index);
        VectorType cutRow();
        VectorType cutColumn();
        DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, Dynamic, maxColumn> cutMatrixRow(size_t from);
        DenseMatrix<T, DenseMatrixType::Row | DenseMatrixType::Vector, maxRow, Dynamic> cutMatrixColumn(size_t from);
        void rowSwap(size_t r1, size_t r2) noexcept;
        void columnSwap(size_t c1, size_t c2) noexcept;
        void rowReduce(size_t r1, size_t r2, size_t element);
        void columnReduce(size_t c1, size_t c2, size_t element);
        /* Getters */
        [[nodiscard]] constexpr static int getType() { return DenseMatrixType::Row | DenseMatrixType::Vector; }
        [[nodiscard]] size_t getRow() const { return Base::getLength(); }
        [[nodiscard]] size_t getColumn() const { return Base::operator[](0).getLength(); }
    };
    /* Operators */
    template<class T, int type, size_t maxRow, size_t maxColumn>
    std::ostream& operator<<(std::ostream& os, const DenseMatrix<T, type, maxRow, maxColumn>& m);

    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> operator+(
            const DenseMatrix<T, type, maxRow, maxColumn>& m1, const DenseMatrix<T, type, maxRow, maxColumn>& m2);

    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> operator-(
            const DenseMatrix<T, type, maxRow, maxColumn>& m1, const DenseMatrix<T, type, maxRow, maxColumn>& m2);

    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> operator*(
            const DenseMatrix<T, type, maxRow, maxColumn>& m1, const DenseMatrix<T, type, maxRow, maxColumn>& m2);

    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> operator*(const DenseMatrix<T, type, maxRow, maxColumn>& m, const MultiScalar& n);
    /* Inline Implementations */
    template<class T, int type, size_t maxRow, size_t maxColumn>
    inline void operator+=(DenseMatrix<T, type, maxRow, maxColumn>& m1
            , const DenseMatrix<T, type, maxRow, maxColumn>& m2) { m1 = m1 + m2; }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    inline void operator-=(DenseMatrix<T, type, maxRow, maxColumn>& m1
            , const DenseMatrix<T, type, maxRow, maxColumn>& m2) { m1 = m1 - m2; }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    inline void operator*=(DenseMatrix<T, type, maxRow, maxColumn>& m
            , const MultiScalar& n) { m = m * n; }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    inline void operator*=(DenseMatrix<T, type, maxRow, maxColumn>& m1
            , const DenseMatrix<T, type, maxRow, maxColumn>& m2) { m1 = m1 * m2; }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    inline void swap(DenseMatrix<T, type, maxRow, maxColumn>& m1
            , DenseMatrix<T, type, maxRow, maxColumn>& m2) noexcept { m1.swap(m2); }
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> reciprocal(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> sqrt(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> factorial(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> ln(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> log(const DenseMatrix<T, type, maxRow, maxColumn>& m, const MultiScalar& a);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> exp(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> pow(const DenseMatrix<T, type, maxRow, maxColumn>& m, const MultiScalar& a);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> cos(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> sin(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> tan(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> sec(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> csc(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> cot(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arccos(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arcsin(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arctan(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arcsec(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arccsc(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arccot(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> cosh(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> sinh(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> tanh(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> sech(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> csch(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> coth(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arccosh(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arcsinh(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arctanh(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arcsech(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arccsch(const DenseMatrix<T, type, maxRow, maxColumn>& m);
    
    template<class T, int type, size_t maxRow, size_t maxColumn>
    DenseMatrix<T, type, maxRow, maxColumn> arccoth(const DenseMatrix<T, type, maxRow, maxColumn>& m);

    namespace Intenal {
        template<class T, int type, size_t maxRow, size_t maxColumn>
        class MatrixTrait<DenseMatrix<T, type, maxRow, maxColumn>> {
        public:
            using ScalarType = T;
            constexpr static size_t MaxRowAtCompile = maxRow;
            constexpr static size_t MaxColumnAtCompile = maxColumn;
        };
    }
}

#include "DenseMatrixImpl.h"

#endif