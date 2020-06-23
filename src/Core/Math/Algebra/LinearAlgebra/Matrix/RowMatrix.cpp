/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "Core/Header/RowMatrix.h"

namespace Physica::Core {
    RowMatrix::RowMatrix() : Matrix(Row) {}

    RowMatrix::RowMatrix(size_t capacity) : Matrix(capacity, Row) {}
    /*!
     * Low level api. Designed for performance.
     * Warning: The first \length elements have not allocated. DO NOT try to visit them.
     */
    RowMatrix::RowMatrix(size_t length, size_t capacity) : Matrix(length, capacity, Row) {}

    RowMatrix::RowMatrix(const CStyleArray<Vector, Dynamic>& array) : Matrix(array, Row) {}

    RowMatrix::RowMatrix(CStyleArray<Vector, Dynamic>&& array) noexcept : Matrix(std::move(array), Row) {}

    RowMatrix::RowMatrix(RowMatrix&& matrix) noexcept : Matrix(std::move(matrix)) {}
    //!Optimize: may be inline append().
    void RowMatrix::appendRow(const Vector& v) {
        Q_ASSERT(v.getLength() == column());
        append(v);
    }

    void RowMatrix::appendRow(Vector&& v) noexcept {
        Q_ASSERT(v.getLength() == column());
        append(std::move(v));
    }

    void RowMatrix::appendColumn(const Vector& v) {
        const auto length = (*this)[0].getLength();
        Q_ASSERT(length == v.getLength());
        for(size_t i = 0; i < length; ++i)
            (*this)[i].append(v[i]);
    }

    void RowMatrix::appendColumn(Vector&& v) noexcept {
        const auto length = (*this)[0].getLength();
        Q_ASSERT(length == v.getLength());
        for(size_t i = 0; i < length; ++i)
            (*this)[i].append(std::move(v[i]));
    }

    void RowMatrix::appendMatrixRow(const Matrix& m) {
        Q_ASSERT((*this)[0].getLength() == m[0].getLength());
        append(m);
    }

    void RowMatrix::appendMatrixRow(Matrix&& m) {
        Q_ASSERT((*this)[0].getLength() == m[0].getLength());
        append(std::move(m));
    }

    void RowMatrix::appendMatrixColumn(Matrix&& m) {
        const auto length = getLength();
        Q_ASSERT(length == m.getLength());
        for(size_t i = 0; i < length; ++i)
            (*this)[i].append(std::move(m[i]));
    }

    void RowMatrix::appendMatrixColumn(const Matrix& m) {
        const auto length = getLength();
        Q_ASSERT(length == m.getLength());
        for(size_t i = 0; i < length; ++i)
            (*this)[i].append(m[i]);
    }

    Vector RowMatrix::cutRow() {
        return cutLast();
    }

    Vector RowMatrix::cutColumn() {
        const auto row = RowMatrix::row();
        Vector result(row, row);
        for(size_t i = 0; i < row; ++i)
            result.allocate((*this)[i].cutLast(), i);
        return result;
    }

    std::unique_ptr<Matrix> RowMatrix::cutMatrixRow(size_t from) {
        return std::unique_ptr<Matrix>(new RowMatrix(cut(from)));
    }

    std::unique_ptr<Matrix> RowMatrix::cutMatrixColumn(size_t from) {
        const auto column = RowMatrix::column();
        auto result = new RowMatrix(column, column);
        for(size_t i = 0; i < column; ++i)
            result->allocate(Vector((*this)[i].cut(from)), i);
        return std::unique_ptr<Matrix>(result);
    }

    void RowMatrix::rowSwap(size_t r1, size_t r2) noexcept {
        Physica::Core::swap((*this)[r1], (*this)[r2]);
    }

    void RowMatrix::columnSwap(size_t c1, size_t c2) noexcept {
        const auto length = getLength();
        for(size_t i = 0; i < length; ++i) {
            auto& row = (*this)[i];
            Physica::Core::swap(row[c1], row[c2]);
        }
    }
    //!Reduce the element at \r2 using \r1
    void RowMatrix::rowReduce(size_t r1, size_t r2, size_t element) {
        Scalar dividend = (*this)(element, r2) / (*this)(element, r1);
        (*this)[r2] -= (*this)[r1] * dividend;
    }
    //!Reduce the element at \c2 using \c1.
    void RowMatrix::columnReduce(size_t c1, size_t c2, size_t element) {
        const auto& element_row = (*this)[element];
        Scalar dividend = element_row[c2] / element_row[c1];
        const auto length = getLength();
        for(size_t i = 0; i < length; ++i) {
            auto& row = (*this)[i];
            row[c2] -= row[c1] * dividend;
        }
    }
}