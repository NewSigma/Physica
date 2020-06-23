/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <iomanip>
#include "Core/Header/Const.h"
#include "Core/Header/Scalar.h"
#include "Core/Header/Matrix.h"
#include "Core/Header/Vector.h"
#include "Core/Header/ColumnMatrix.h"
#include "Core/Header/RowMatrix.h"

namespace Physica::Core {
    ////////////////////////////////////////Column Matrix////////////////////////////////////////////
    Matrix::Matrix(MatrixType type) : CStyleArray<Vector, Dynamic>(), type(type) {}

    Matrix::Matrix(size_t capacity, MatrixType type) : CStyleArray<Vector, Dynamic>(capacity), type(type) {}
    /*!
     * Low level api. Designed for performance.
     * Warning: The first \length elements have not allocated. DO NOT try to visit them.
     */
    Matrix::Matrix(size_t length, size_t capacity, MatrixType type)
            : CStyleArray<Vector, Dynamic>(length, capacity), type(type) {}

    Matrix::Matrix(const CStyleArray<Vector, Dynamic>& array, MatrixType type)
            : CStyleArray<Vector, Dynamic>(array), type(type) {}

    Matrix::Matrix(CStyleArray<Vector, Dynamic>&& array, MatrixType type) noexcept
            : CStyleArray<Vector, Dynamic>(std::move(array)), type(type) {}

    Matrix::Matrix(Matrix&& matrix) noexcept : CStyleArray<Vector, Dynamic>(std::move(matrix)), type(matrix.type) {}

    Matrix& Matrix::operator=(const Matrix& m) {
        if(this == &m)
            return *this;
        this->~Matrix();
        type = m.type;
        CStyleArray<Vector, Dynamic>::operator=(m);
        return *this;
    }

    Matrix& Matrix::operator=(Matrix&& m) noexcept {
        this->~Matrix();
        type = m.type;
        CStyleArray<Vector, Dynamic>::operator=(std::move(m));
        return *this;
    }
    /*!
     * Eliminate elements at column index using row reduce, which are above the row index.
     * Warning: Element at (index, index) must not be zero. Or a divide by zero exception will be thrown.
     */
    void Matrix::upperEliminate(size_t index) {
        for(size_t i = 0; i < index; ++i)
            rowReduce(index, i, index);
    }
    /*!
     * Eliminate elements at column index using row reduce, which are above the row index.
     * Warning: Element at (index, index) must not be zero. Or a divide by zero exception will be thrown.
     */
    void Matrix::lowerEliminate(size_t index) {
        const auto r = row();
        for(size_t i = index + 1; i < r; ++i)
            rowReduce(index, i, index);
    }

    void Matrix::impactPivoting() {
        const auto r = row();
        for(size_t i = 0; i < r; ++i)
            impactPivoting(i);
    }

    void Matrix::swap(Matrix& m) noexcept {
        CStyleArray<Vector, Dynamic>::swap(m);
        auto temp = type;
        type = m.type;
        m.type = temp;
    }
    //Print all elements.
    std::ostream& operator<<(std::ostream& os, const Matrix& m) {
        const auto row = m.row();
        const auto column = m.column();
        //10 is the max precision of double.
        os << std::setprecision(10);
        for(size_t i = 0; i < row; ++i) {
            for(size_t j = 0; j < column; ++j)
                os << std::to_string(double(m(i, j))) << '\t';
            os << '\n';
        }
        //6 is the default precision.
        os << std::setprecision(6);
        return os;
    }

    std::unique_ptr<Matrix> operator+(const Matrix& m1, const Matrix& m2) {
        Q_ASSERT(m1.row() == m2.row() && m1.column() == m2.column());
        const auto length = m1.getLength();
        std::unique_ptr<Matrix> result(
                m1.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(m1[i] + m2[i], i);
        return result;
    }

    std::unique_ptr<Matrix> operator-(const Matrix& m1, const Matrix& m2) {
        Q_ASSERT(m1.row() == m2.row() && m1.column() == m2.column());
        const auto length = m1.getLength();
        std::unique_ptr<Matrix> result(
                m1.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(m1[i] - m2[i], i);
        return result;
    }

    std::unique_ptr<Matrix> operator*(const Matrix& m, const Scalar& n) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(m[i] * n, i);
        return result;
    }

    std::unique_ptr<Matrix> operator*(const Matrix& m1, const Matrix& m2) {
        Q_ASSERT(m1.column() == m2.row());
        const bool isColumn = m1.getType() == Matrix::Column;
        const auto result_row = m1.row();
        const auto result_column = m2.column();
        const auto m1_column = m1.column();
        const auto matrix_size = isColumn ? result_column : result_row;
        const auto vector_length = isColumn ? result_row : result_column;
        std::unique_ptr<Matrix> result(
                isColumn
                ? static_cast<Matrix*>(new ColumnMatrix(matrix_size, matrix_size))
                : static_cast<Matrix*>(new RowMatrix(matrix_size, matrix_size)));

        size_t i, j;
        const size_t& r = isColumn ? j : i;
        const size_t& c = isColumn ? i : j;
        for(i = 0; i < matrix_size; ++i) {
            Vector new_vector(vector_length, vector_length);
            for(j = 0; j < vector_length; ++j) {
                Scalar element(BasicConst::getInstance().get_0());
                for(size_t k = 0; k < m1_column; ++k)
                    element += m1(r, k) * m2(k, c);
                new_vector.allocate(std::move(element), j);
            }
        }
        return result;
    }
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    std::unique_ptr<Matrix> reciprocal(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(reciprocal(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> sqrt(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(sqrt(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> factorial(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(factorial(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> ln(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(ln(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> log(const Matrix& m, const Scalar& a) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(log(m[i], a), i);
        return result;
    }

    std::unique_ptr<Matrix> exp(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(exp(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> pow(const Matrix& m, const Scalar& a) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(pow(m[i], a), i);
        return result;
    }

    std::unique_ptr<Matrix> cos(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(cos(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> sin(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(sin(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> tan(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(tan(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> sec(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(sec(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> csc(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(csc(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> cot(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(cot(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> arccos(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(arccos(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> arcsin(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(arcsin(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> arctan(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(arctan(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> arcsec(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(arcsec(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> arccsc(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(arccsc(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> arccot(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(arccot(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> cosh(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(cosh(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> sinh(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(sinh(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> tanh(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(tanh(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> sech(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(sech(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> csch(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(csch(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> coth(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(coth(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> arccosh(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(arccosh(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> arcsinh(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(arcsinh(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> arctanh(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(arctanh(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> arcsech(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(arcsech(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> arccsch(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(arccsch(m[i]), i);
        return result;
    }

    std::unique_ptr<Matrix> arccoth(const Matrix& m) {
        const auto length = m.getLength();
        std::unique_ptr<Matrix> result(
                m.getType() == Matrix::Column
                ? static_cast<Matrix*>(new ColumnMatrix(length, length))
                : static_cast<Matrix*>(new RowMatrix(length, length)));
        for(size_t i = 0; i < length; ++i)
            (*result).allocate(arccoth(m[i]), i);
        return result;
    }
}