/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Core/Header/Const.h"
#include "Core/Header/Numerical.h"
#include "Core/Header/Matrix.h"
#include "Core/Header/Vector.h"
#include <iomanip>

namespace Physica::Core {
    ////////////////////////////////////////Column Matrix////////////////////////////////////////////
    Matrix::Matrix(MatrixType type) : vectors(nullptr), length(0), capacity(0), type(type) {}

    Matrix::Matrix(Vector* arr, size_t length, MatrixType type)
            : vectors(arr), length(length), capacity(length), type(type) {
        arr = nullptr;
    }

    Matrix::Matrix(const Matrix& matrix)
            : vectors(new Vector[matrix.length]), length(matrix.length), capacity(matrix.capacity), type(matrix.type) {
        for(size_t i = 0; i < length; ++i)
            vectors[i] = matrix.vectors[i];
    }

    Matrix::Matrix(Matrix&& matrix) noexcept
            : vectors(matrix.vectors), length(matrix.length), capacity(matrix.capacity), type(matrix.type) {
        matrix.vectors = nullptr;
    }

    Matrix::~Matrix() {
        delete[] vectors;
    }

    Matrix& Matrix::operator=(const Matrix& m) noexcept {
        if(this == &m)
            return *this;
        this->~Matrix();
        vectors = new Vector[m.length];
        capacity = m.capacity;
        length = m.length;
        for(size_t i = 0; i < length; ++i)
            vectors[i] = m.vectors[i];
        return *this;
    }

    Matrix& Matrix::operator=(Matrix&& m) noexcept {
        this->~Matrix();
        vectors = m.vectors;
        length = m.length;
        capacity = m.capacity;
        m.vectors = nullptr;
        return *this;
    }
    /*
     * Eliminate elements at column index using row reduce, which are above the row index.
     * Warning: Element at (index, index) must not be zero. Or a divide by zero exception will be thrown.
     */
    void Matrix::upperEliminate(size_t index) {
        for(size_t i = 0; i < index; ++i)
            rowReduce(index, i, index);
    }
    /*
     * Eliminate elements at column index using row reduce, which are above the row index.
     * Warning: Element at (index, index) must not be zero. Or a divide by zero exception will be thrown.
     */
    void Matrix::lowerEliminate(size_t index) {
        const auto col = column();
        for(size_t i = index + 1; i < col; ++i)
            rowReduce(index, i, index);
    }

    void Matrix::impactPivoting() {
        const auto r = row();
        for(size_t i = 0; i < r; ++i)
            impactPivoting(i);
    }

    void Matrix::swap(Matrix& m) noexcept {
        auto temp = vectors;
        vectors = m.vectors;
        m.vectors = temp;
        auto temp_size = length;
        length = m.length;
        m.length = temp_size;
        temp_size = capacity;
        capacity = m.capacity;
        m.capacity = temp_size;
    }
    //Print all vectors.
    std::ostream& operator<<(std::ostream& os, const Matrix& m) {
        const auto row = m.row();
        const auto column = m.column();
        //10 is the max precision of double.
        os << std::setprecision(10);
        for(size_t i = 0; i < column; ++i) {
            for(size_t j = 0; j < row; ++j)
                os << std::to_string(double(m(j, i))) << '\t';
            os << '\n';
        }
        //6 is the default precision.
        os << std::setprecision(6);
        return os;
    }

    std::unique_ptr<Matrix> operator+(const Matrix& m1, const Matrix& m2) {
        Q_ASSERT(m1.row() == m2.row());
        Q_ASSERT(m1.column() == m2.column());
        const auto length = m1.getLength();
        auto new_vectors = new Vector[length];
        for(size_t i = 0; i < length; ++i)
            new_vectors[i] = m1[i] + m2[i];
        return std::unique_ptr<Matrix>(
                m1.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(new_vectors, length))
                : static_cast<Matrix*>(new RowMatrix(new_vectors, length)));
    }

    std::unique_ptr<Matrix> operator-(const Matrix& m1, const Matrix& m2) {
        Q_ASSERT(m1.row() == m2.row());
        Q_ASSERT(m1.column() == m2.column());
        const auto length = m1.getLength();
        auto new_vectors = new Vector[length];
        for(size_t i = 0; i < length; ++i)
            new_vectors[i] = m1[i] - m2[i];
        return std::unique_ptr<Matrix>(
                m1.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(new_vectors, length))
                : static_cast<Matrix*>(new RowMatrix(new_vectors, length)));
    }

    std::unique_ptr<Matrix> operator*(const Matrix& m, const Numerical& n) {
        const auto length = m.getLength();
        auto new_vectors = new Vector[length];
        for(size_t i = 0; i < length; ++i)
            new_vectors[i] = m[i] * n;
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(new_vectors, length))
                : static_cast<Matrix*>(new RowMatrix(new_vectors, length)));
    }

    std::unique_ptr<Matrix> operator*(const Matrix& m1, const Matrix& m2) {
        Q_ASSERT(m1.column() == m2.row());
        const auto result_row = m1.row();
        const auto result_column = m2.column();
        const auto m1_column = m1.column();
        auto new_vectors = new Vector[result_column];
        for(size_t i = 0; i < result_column; ++i)
            new_vectors[i].initVector(result_row);
        for(size_t i = 0; i < result_column; ++i) {
            for(size_t j = 0; j < result_row; ++j) {
                auto& element = new_vectors[i][j];
                element = BasicConst::getInstance().get_0();
                for(size_t k = 0; k < m1_column; ++k)
                    element += m1(j, k) * m2(k, i);
            }
        }
        return std::unique_ptr<Matrix>(
                m1.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(new_vectors, result_column))
                : static_cast<Matrix*>(new RowMatrix(new_vectors, result_row)));
    }
    ////////////////////////////////////////Column Matrix////////////////////////////////////////////
    ColumnMatrix::ColumnMatrix() : Matrix(Column) {}

    ColumnMatrix::ColumnMatrix(size_t column, size_t row) : Matrix(new Vector[column], column, Column) {
        for(size_t i = 0; i < column; ++i)
            vectors[i].initVector(row);
    }

    ColumnMatrix::ColumnMatrix(Vector* vectors, size_t length) : Matrix(vectors, length, Column) {}
    //Warning: Length of v must equals to column() or a out of bounder visit will happen.
    void ColumnMatrix::appendRow(Vector v) {
        indirectVectorAppend(std::move(v));
    }
    //Warning: Length of v must equals to row() or a out of bounder visit will happen.
    void ColumnMatrix::appendColumn(Vector v) {
        directVectorAppend(std::move(v));
    }
    //Warning: m must be a ColumnMatrix and both of the two length must be equal.
    void ColumnMatrix::appendMatrixRow(Matrix&& m) {
        indirectMatrixAppend(std::move(m));
    }
    //Warning: m must be a ColumnMatrix and both of the two length must be equal.
    void ColumnMatrix::appendMatrixColumn(Matrix&& m) {
        directMatrixAppend(std::move(m));
    }

    Vector ColumnMatrix::cutRow() {
        return indirectVectorCut();
    }

    Vector ColumnMatrix::cutColumn() {
        return directVectorCut();
    }

    std::unique_ptr<Matrix> ColumnMatrix::cutMatrixRow(size_t from) {
        return std::unique_ptr<Matrix>(indirectMatrixCut(from));
    }

    std::unique_ptr<Matrix> ColumnMatrix::cutMatrixColumn(size_t from) {
        return std::unique_ptr<Matrix>(directMatrixCut(from));
    }

    void ColumnMatrix::rowSwap(size_t r1, size_t r2) noexcept {
        indirectSwap(r1, r2);
    }

    void ColumnMatrix::columnSwap(size_t c1, size_t c2) noexcept {
        directSwap(c1, c2);
    }

    void ColumnMatrix::rowReduce(size_t r1, size_t r2, size_t element) {
        const auto r = row();
        Q_ASSERT(r1 < r && r2 < r && element < r);
        indirectReduce(r1, r2, element);
    }

    void ColumnMatrix::columnReduce(size_t c1, size_t c2, size_t element) {
        const auto c = column();
        Q_ASSERT(c1 < c && c2 < c && element < c);
        directReduce(c1, c2, element);
    }
    ////////////////////////////////////////Row Matrix////////////////////////////////////////////
    RowMatrix::RowMatrix() : Matrix(Row) {}

    RowMatrix::RowMatrix(size_t column, size_t row) : Matrix(new Vector[row], row, Row) {
        for(size_t i = 0; i < row; ++i)
            vectors[i].initVector(column);
    }

    RowMatrix::RowMatrix(Vector* vectors, size_t length) : Matrix(vectors, length, Row) {}
    //Warning: Length of v must equals to column() or a out of bounder visit will happen.
    void RowMatrix::appendRow(Vector v) {
        directVectorAppend(std::move(v));
    }
    //Warning: Length of v must equals to row() or a out of bounder visit will happen.
    void RowMatrix::appendColumn(Vector v) {
        indirectVectorAppend(std::move(v));
    }
    //Warning: m must be a RowMatrix and both of the two length must be equal.
    void RowMatrix::appendMatrixRow(Matrix&& m) {
        directMatrixAppend(std::move(m));
    }
    //Warning: m must be a RowMatrix and both of the two length must be equal.
    void RowMatrix::appendMatrixColumn(Matrix&& m) {
        indirectMatrixAppend(std::move(m));
    }

    Vector RowMatrix::cutRow() {
        return directVectorCut();
    }

    Vector RowMatrix::cutColumn() {
        return indirectVectorCut();
    }

    std::unique_ptr<Matrix> RowMatrix::cutMatrixRow(size_t from) {
        return std::unique_ptr<Matrix>(directMatrixCut(from));
    }

    std::unique_ptr<Matrix> RowMatrix::cutMatrixColumn(size_t from) {
        return std::unique_ptr<Matrix>(indirectMatrixCut(from));
    }

    void RowMatrix::rowSwap(size_t r1, size_t r2) noexcept {
        directSwap(r1, r2);
    }

    void RowMatrix::columnSwap(size_t c1, size_t c2) noexcept {
        indirectSwap(c1, c2);
    }

    void RowMatrix::rowReduce(size_t r1, size_t r2, size_t element) {
        const auto r = row();
        Q_ASSERT(r1 < r && r2 < r && element < r);
        directReduce(r1, r2, element);
    }

    void RowMatrix::columnReduce(size_t c1, size_t c2, size_t element) {
        const auto c = column();
        Q_ASSERT(c1 < c && c2 < c && element < c);
        indirectReduce(c1, c2, element);
    }
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    std::unique_ptr<Matrix> reciprocal(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = reciprocal(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> sqrt(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = sqrt(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> factorial(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = factorial(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> ln(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = ln(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> log(const Matrix& m, const Numerical& a) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = log(m[i], a);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> exp(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = exp(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> pow(const Matrix& m, const Numerical& a) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = pow(m[i], a);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> cos(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = cos(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> sin(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = sin(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> tan(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = tan(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> sec(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = sec(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> csc(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = csc(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> cot(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = cot(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arccos(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arccos(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arcsin(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arcsin(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arctan(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arctan(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arcsec(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arcsec(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arccsc(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arccsc(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arccot(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arccot(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> cosh(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = cosh(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> sinh(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = sinh(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> tanh(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = tanh(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> sech(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = sech(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> csch(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = csch(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> coth(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = coth(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arccosh(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arccosh(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arcsinh(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arcsinh(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arctanh(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arctanh(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arcsech(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arcsech(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arccsch(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arccsch(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arccoth(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arccoth(m[i]);
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }
}