/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Core/Header/Const.h"
#include "Core/Header/Numerical.h"
#include "Core/Header/Matrix.h"
#include "Core/Header/Vector.h"

namespace Physica::Core {
    Matrix::Matrix() : vectors(nullptr), length(0) {}

    Matrix::Matrix(size_t c, size_t r) : Matrix(new Vector[c], c) {
        for(size_t i = 0; i < c; ++i)
            vectors[i].initVector(r);
    }

    Matrix::Matrix(Vector*& arr, size_t length) : vectors(arr), length(length) {
        arr = nullptr;
    }

    Matrix::Matrix(Vector*&& arr, size_t length) : vectors(arr), length(length) {}

    Matrix::Matrix(Matrix& matrix) : vectors(new Vector[matrix.length]), length(matrix.length)  {
        for(size_t i = 0; i < length; ++i)
            vectors[i] = matrix.vectors[i];
    }

    Matrix::Matrix(Matrix&& matrix) noexcept : vectors(matrix.vectors), length(matrix.length) {
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
        length = m.length;
        for(size_t i = 0; i < length; ++i)
            vectors[i] = m.vectors[i];
        return *this;
    }

    Matrix& Matrix::operator=(Matrix&& m) noexcept {
        this->~Matrix();
        vectors = m.vectors;
        length = m.length;
        m.vectors = nullptr;
        return *this;
    }

    Matrix operator+(const Matrix& m1, const Matrix& m2) {
        auto new_vectors = new Vector[m1.getLength()];
        for(size_t i = 0; i < m1.getLength(); ++i)
            new_vectors[i] = m1[i] + m2[i];
        return Matrix(new_vectors, m1.getLength());
    }

    Matrix operator-(const Matrix& m1, const Matrix& m2) {
        auto new_vectors = new Vector[m1.getLength()];
        for(size_t i = 0; i < m1.getLength(); ++i)
            new_vectors[i] = m1[i] - m2[i];
        return Matrix(new_vectors, m1.getLength());
    }

    Matrix operator*(const Matrix& m1, const Matrix& m2) {
        const auto result_row = m1.row();
        const auto result_column = m2.column();
        auto new_vectors = new Vector[result_column];
        for(size_t i = 0; i < result_column; ++i)
            new_vectors[i].initVector(result_row);
        for(size_t i = 0; i < result_column; ++i) {
            for(size_t j = 0; j < result_row; ++j) {
                auto& element = new_vectors[i][j];
                element = BasicConst::getInstance().get_0();
                for(size_t k = 0; k < result_column; ++k)
                    element += m1(j, k) * m2(k, i);
            }
        }
        return Matrix(new_vectors, m1.getLength());
    }

    Matrix operator*(const Matrix& m, const Numerical& n) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = m[i] * n;
        return Matrix(arr, m.getLength());
    }
    //Warning: Indexes start from 0, should avoid out of border access.
    void Matrix::rowSwap(size_t r1, size_t r2) noexcept {
        for(size_t i = 0; i < length; ++i)
            swap(vectors[i][r1], vectors[i][r2]);
    }
    //Warning: Indexes start from 0, should avoid out of border access.
    void Matrix::columnSwap(size_t c1, size_t c2) noexcept {
        auto temp = vectors[c1];
        vectors[c1] = vectors[c2];
        vectors[c2] = temp;
    }
    /*
     * Eliminate the element at r2 using r1.
     * Warning:
     * 1.Indexes start from 0, should avoid out of border access.
     * 2.vectors[element][r1] mustn't be zero.
     */
    void Matrix::rowEliminate(size_t r1, size_t r2, size_t element) {
        auto& elementColumnVector = vectors[element];
        Numerical dividend = elementColumnVector[r2] / elementColumnVector[r1];
        for(size_t i = 0; i < length; ++i)
            vectors[i][r2] -= vectors[i][r1] * dividend;
    }
    /*
     * Eliminate the element at c2 using c1.
     * Warning:
     * 1.Indexes start from 0, should avoid out of border access.
     * 2.vectors[c1][element] mustn't be zero.
     */
    void Matrix::columnEliminate(size_t c1, size_t c2, size_t element) {
        auto& columnVector1 = vectors[c1];
        auto& columnVector2 = vectors[c2];
        Numerical dividend = columnVector2[element] / columnVector1[element];
        const auto matrixRow = row();
        for(size_t i = 0; i < matrixRow; ++i)
            columnVector2[i] -= columnVector1[i] * dividend;
    }

    void Matrix::transpose() {
        const size_t new_row = column();
        const size_t new_column = row();
        auto arr = new Vector[new_column];
        for(size_t i = 0; i < new_column; ++i)
            arr[i].initVector(new_row);
        for(size_t i = 0; i < new_row; ++i) {
            for(size_t j = 0; j < new_column; ++j)
                arr[j][i] = std::move(vectors[i][j]);
        }
        this->~Matrix();
        length = new_column;
        vectors = arr;
    }
////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    Matrix reciprocal(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = reciprocal(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix sqrt(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = sqrt(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix factorial(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = factorial(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix ln(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = ln(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix log(const Matrix& m, const Numerical& a) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = log(m[i], a);
        return Matrix(arr, m.getLength());
    }

    Matrix exp(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = exp(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix pow(const Matrix& m, const Numerical& a) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = pow(m[i], a);
        return Matrix(arr, m.getLength());
    }

    Matrix cos(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = cos(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix sin(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = sin(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix tan(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = tan(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix sec(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = sec(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix csc(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = csc(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix cot(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = cot(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arccos(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arccos(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arcsin(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arcsin(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arctan(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arctan(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arcsec(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arcsec(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arccsc(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arccsc(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arccot(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arccot(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix cosh(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = cosh(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix sinh(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = sinh(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix tanh(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = tanh(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix sech(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = sech(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix csch(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = csch(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix coth(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = coth(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arccosh(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arccosh(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arcsinh(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arcsinh(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arctanh(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arctanh(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arcsech(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arcsech(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arccsch(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arccsch(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arccoth(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(size_t i = 0; i < m.getLength(); ++i)
            arr[i] = arccoth(m[i]);
        return Matrix(arr, m.getLength());
    }
}