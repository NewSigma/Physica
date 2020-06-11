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

    Matrix::Matrix(size_t length, MatrixType type)
            : vectors(reinterpret_cast<Vector*>(malloc(length * sizeof(Vector))))
            , length(0), capacity(length), type(type) {}
    //Declared as protected to avoid incorrect use.
    Matrix::Matrix(Vector* arr, size_t length, MatrixType type)
            : vectors(arr), length(length), capacity(length), type(type) {}

    Matrix::Matrix(const Matrix& matrix)
            : vectors(reinterpret_cast<Vector*>(malloc(matrix.capacity * sizeof(Vector))))
            , length(matrix.length), capacity(matrix.capacity), type(matrix.type) {
        for(size_t i = 0; i < length; ++i)
            new (vectors + i) Vector(matrix.vectors[i]);
    }

    Matrix::Matrix(Matrix&& matrix) noexcept
            : vectors(matrix.vectors), length(matrix.length), capacity(matrix.capacity), type(matrix.type) {
        matrix.vectors = nullptr;
        matrix.length = 0;
    }

    Matrix::~Matrix() {
        for(size_t i = 0; i < length; ++i)
            (vectors + i)->~Vector();
        free(vectors);
    }

    Matrix& Matrix::operator=(const Matrix& m) noexcept {
        if(this == &m)
            return *this;
        this->~Matrix();
        length = m.length;
        capacity = m.capacity;
        vectors = reinterpret_cast<Vector*>(malloc(capacity * sizeof(Vector)));
        for(size_t i = 0; i < length; ++i)
            new (vectors + i) Vector(m.vectors[i]);
        return *this;
    }

    Matrix& Matrix::operator=(Matrix&& m) noexcept {
        this->~Matrix();
        vectors = m.vectors;
        length = m.length;
        capacity = m.capacity;
        m.vectors = nullptr;
        m.length = 0;
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

    void Matrix::resize(size_t size) {
        if(length > size) {
            for(size_t i = size; size < length; ++i)
                (vectors + i)->~Vector();
            length = size;
        }
        vectors = reinterpret_cast<Vector*>(realloc(vectors, size * sizeof(Vector)));
        capacity = size;
    }

    void Matrix::squeeze() {
        vectors = reinterpret_cast<Vector*>(realloc(vectors, length * sizeof(Vector)));
        capacity = length;
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
        Q_ASSERT(m1.row() == m2.row());
        Q_ASSERT(m1.column() == m2.column());
        const auto length = m1.getLength();
        auto new_vectors = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < length; ++i)
            new (new_vectors + i) Vector(m1[i] + m2[i]);
        return std::unique_ptr<Matrix>(
                m1.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(new_vectors, length))
                : static_cast<Matrix*>(new RowMatrix(new_vectors, length)));
    }

    std::unique_ptr<Matrix> operator-(const Matrix& m1, const Matrix& m2) {
        Q_ASSERT(m1.row() == m2.row());
        Q_ASSERT(m1.column() == m2.column());
        const auto length = m1.getLength();
        auto new_vectors = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < length; ++i)
            new (new_vectors + i) Vector(m1[i] - m2[i]);
        return std::unique_ptr<Matrix>(
                m1.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(new_vectors, length))
                : static_cast<Matrix*>(new RowMatrix(new_vectors, length)));
    }

    std::unique_ptr<Matrix> operator*(const Matrix& m, const Numerical& n) {
        const auto length = m.getLength();
        auto new_vectors = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < length; ++i)
            new (new_vectors + i) Vector(m[i] * n);
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
        auto new_vectors = reinterpret_cast<Vector*>(malloc(result_column * sizeof(Vector)));
        for(size_t i = 0; i < result_column; ++i) {
            auto new_vector = new (new_vectors + i) Vector();
            for(size_t j = 0; j < result_row; ++j) {
                Numerical element(BasicConst::getInstance().get_0());
                for(size_t k = 0; k < m1_column; ++k)
                    element += m1(j, k) * m2(k, i);
                new_vector->grow(std::move(element));
            }
        }
        return std::unique_ptr<Matrix>(
                m1.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(new_vectors, result_column))
                : static_cast<Matrix*>(new RowMatrix(new_vectors, result_row)));
    }
    ////////////////////////////////////////Column Matrix////////////////////////////////////////////
    ColumnMatrix::ColumnMatrix() : Matrix(Column) {}

    ColumnMatrix::ColumnMatrix(size_t column, size_t row)
            : Matrix(reinterpret_cast<Vector*>(malloc(column * sizeof(Vector))), column, Column) {
        for(size_t i = 0; i < column; ++i)
            new (vectors + i) Vector(row);
    }

    ColumnMatrix::ColumnMatrix(Vector* vectors, size_t length) : Matrix(vectors, length, Column) {}

    ColumnMatrix::ColumnMatrix(ColumnMatrix&& matrix) noexcept : Matrix(std::move(matrix)) {}
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
        indirectReduce(r1, r2, element);
    }

    void ColumnMatrix::columnReduce(size_t c1, size_t c2, size_t element) {
        directReduce(c1, c2, element);
    }
    ////////////////////////////////////////Row Matrix////////////////////////////////////////////
    RowMatrix::RowMatrix() : Matrix(Row) {}

    RowMatrix::RowMatrix(size_t column, size_t row)
            : Matrix(reinterpret_cast<Vector*>(malloc(row * sizeof(Vector))), row, Row) {
        for(size_t i = 0; i < row; ++i)
            new (vectors + i) Vector(column);
    }

    RowMatrix::RowMatrix(Vector* vectors, size_t length) : Matrix(vectors, length, Row) {}

    RowMatrix::RowMatrix(RowMatrix&& matrix) noexcept : Matrix(std::move(matrix)) {}
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
        directReduce(r1, r2, element);
    }

    void RowMatrix::columnReduce(size_t c1, size_t c2, size_t element) {
        indirectReduce(c1, c2, element);
    }
    ////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    std::unique_ptr<Matrix> reciprocal(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(reciprocal(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> sqrt(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(sqrt(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> factorial(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(factorial(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> ln(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(ln(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> log(const Matrix& m, const Numerical& a) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(log(m[i], a));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> exp(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(exp(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> pow(const Matrix& m, const Numerical& a) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(pow(m[i], a));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> cos(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(cos(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> sin(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(sin(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> tan(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(tan(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> sec(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(sec(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> csc(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(csc(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> cot(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(cot(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arccos(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(arccos(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arcsin(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(arcsin(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arctan(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(arctan(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arcsec(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(arcsec(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arccsc(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(arccsc(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arccot(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(arccot(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> cosh(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(cosh(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> sinh(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(sinh(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> tanh(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(tanh(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> sech(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(sech(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> csch(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(csch(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> coth(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(coth(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arccosh(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(arccosh(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arcsinh(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(arcsinh(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arctanh(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(arctanh(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arcsech(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(arcsech(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arccsch(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(arccsch(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }

    std::unique_ptr<Matrix> arccoth(const Matrix& m) {
        const auto length = m.getLength();
        auto arr = reinterpret_cast<Vector*>(malloc(length * sizeof(Vector)));
        for(size_t i = 0; i < m.getLength(); ++i)
            new (arr + i) Vector(arccoth(m[i]));
        return std::unique_ptr<Matrix>(
                m.getType() == Matrix::Column ?
                static_cast<Matrix*>(new ColumnMatrix(arr, m.getLength()))
                : static_cast<Matrix*>(new RowMatrix(arr, m.getLength())));
    }
}