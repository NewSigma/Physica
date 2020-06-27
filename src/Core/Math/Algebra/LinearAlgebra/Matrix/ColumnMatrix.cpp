/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/ColumnMatrix.h"

namespace Physica::Core {
    ColumnMatrix::ColumnMatrix() : Matrix(Column) {}

    ColumnMatrix::ColumnMatrix(size_t capacity) : Matrix(capacity, Column) {}
    /*!
     * Low level api. Designed for performance.
     * Warning: The first \length elements have not allocated. DO NOT try to visit them.
     */
    ColumnMatrix::ColumnMatrix(size_t length, size_t capacity) : Matrix(length, capacity, Column) {}

    ColumnMatrix::ColumnMatrix(const CStyleArray<Vector, Dynamic>& array) : Matrix(array, Column) {}

    ColumnMatrix::ColumnMatrix(CStyleArray<Vector, Dynamic>&& array) noexcept : Matrix(std::move(array), Column) {}

    ColumnMatrix::ColumnMatrix(ColumnMatrix&& matrix) noexcept : Matrix(std::move(matrix)) {}
    //!Optimize: may be inline append().
    void ColumnMatrix::appendRow(const Vector& v) {
        const auto length = (*this)[0].getLength();
        Q_ASSERT(length == v.getLength());
        for(size_t i = 0; i < length; ++i)
            (*this)[i].append(v[i]);
    }

    void ColumnMatrix::appendRow(Vector&& v) noexcept {
        const auto length = (*this)[0].getLength();
        Q_ASSERT(length == v.getLength());
        for(size_t i = 0; i < length; ++i)
            (*this)[i].append(std::move(v[i]));
    }

    void ColumnMatrix::appendColumn(const Vector& v) {
        Q_ASSERT(v.getLength() == row());
        append(v);
    }

    void ColumnMatrix::appendColumn(Vector&& v) noexcept {
        Q_ASSERT(v.getLength() == row());
        append(std::move(v));
    }

    void ColumnMatrix::appendMatrixRow(const Matrix& m) {
        const auto length = getLength();
        Q_ASSERT(length == m.getLength());
        for(size_t i = 0; i < length; ++i)
            (*this)[i].append(m[i]);
    }

    void ColumnMatrix::appendMatrixRow(Matrix&& m) {
        const auto length = getLength();
        Q_ASSERT(length == m.getLength());
        for(size_t i = 0; i < length; ++i)
            (*this)[i].append(std::move(m[i]));
    }

    void ColumnMatrix::appendMatrixColumn(const Matrix& m) {
        Q_ASSERT((*this)[0].getLength() == m[0].getLength());
        append(m);
    }

    void ColumnMatrix::appendMatrixColumn(Matrix&& m) {
        Q_ASSERT((*this)[0].getLength() == m[0].getLength());
        append(std::move(m));
    }

    Vector ColumnMatrix::cutRow() {
        const auto column = ColumnMatrix::column();
        Vector result(column, column);
        for(size_t i = 0; i < column; ++i)
            result.allocate((*this)[i].cutLast(), i);
        return result;
    }

    Vector ColumnMatrix::cutColumn() {
        return cutLast();
    }

    std::unique_ptr<Matrix> ColumnMatrix::cutMatrixRow(size_t from) {
        const auto row = ColumnMatrix::row();
        auto result = new ColumnMatrix(row, row);
        for(size_t i = 0; i < row; ++i)
            result->allocate(Vector((*this)[i].cut(from)), i);
        return std::unique_ptr<Matrix>(result);
    }

    std::unique_ptr<Matrix> ColumnMatrix::cutMatrixColumn(size_t from) {
        return std::unique_ptr<Matrix>(new ColumnMatrix(cut(from)));
    }

    void ColumnMatrix::rowSwap(size_t r1, size_t r2) noexcept {
        const auto length = getLength();
        for(size_t i = 0; i < length; ++i) {
            auto& column = (*this)[i];
            Physica::Core::swap(column[r1], column[r2]);
        }
    }

    void ColumnMatrix::columnSwap(size_t c1, size_t c2) noexcept {
        Physica::Core::swap((*this)[c1], (*this)[c2]);
    }
    //!Reduce the element at \r2 using \r1.
    void ColumnMatrix::rowReduce(size_t r1, size_t r2, size_t element) {
        const auto& element_column = (*this)[element];
        Scalar dividend = element_column[r2] / element_column[r1];
        const auto length = getLength();
        for(size_t i = 0; i < length; ++i) {
            auto& column = (*this)[i];
            column[r2] -= column[r1] * dividend;
        }
    }
    //!Reduce the element at \c2 using \c1.
    void ColumnMatrix::columnReduce(size_t c1, size_t c2, size_t element) {
        Scalar dividend = (*this)(element, c2) / (*this)(element, c1);
        (*this)[c2] -= (*this)[c1] * dividend;
    }
}