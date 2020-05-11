/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Core/Header/Const.h"
#include "Core/Header/Numerical.h"
#include "Core/Header/Matrix.h"
#include "Core/Header/Vector.h"

namespace Physica::Core {
    Matrix::Matrix() : vectors(nullptr), length(0) {}

    Matrix::Matrix(int c, int r) : Matrix(new Vector[c], c) {
        for(int i = 0; i < c; ++i)
            vectors[i] = Vector(r);
    }

    Matrix::Matrix(Vector*& arr, int length) : vectors(arr), length(length) {
        arr = nullptr;
    }

    Matrix::Matrix(Vector*&& arr, int length) : vectors(arr), length(length) {}

    Matrix::Matrix(Matrix& matrix) : vectors(new Vector[matrix.length]), length(matrix.length)  {
        for(int i = 0; i < length; ++i)
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
        for(int i = 0; i < length; ++i)
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
        for(int i = 0; i < m1.getLength(); ++i)
            new_vectors[i] = m1[i] + m2[i];
        return Matrix(new_vectors, m1.getLength());
    }

    Matrix operator-(const Matrix& m1, const Matrix& m2) {
        auto new_vectors = new Vector[m1.getLength()];
        for(int i = 0; i < m1.getLength(); ++i)
            new_vectors[i] = m1[i] - m2[i];
        return Matrix(new_vectors, m1.getLength());
    }

    Matrix operator*(const Matrix& m1, const Matrix& m2) {
        auto new_vectors = new Vector[m2.getLength()];
        for(int i = 0; i < m1.getLength(); ++i)
            new_vectors[i] = Vector(m1.row());
        for(int i = 0; i < m2.getLength(); ++i) {
            for(int j = 0; j < m1.row(); ++j) {
                auto& element = new_vectors[i][j];
                element = BasicConst::getInstance().get_0();
                for(int k = 0; k < m1.getLength(); ++k)
                    element += m1(j, k) * m2(k, i);
            }
        }
        return Matrix(new_vectors, m1.getLength());
    }

    Matrix operator*(const Matrix& m, const Numerical& n) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = m[i] * n;
        return Matrix(arr, m.getLength());
    }

    void Matrix::transpose() {
        auto arr = new Vector[row()];
        for(int i = 0; i < row(); ++i)
            arr[i] = Vector(getLength());
        for(int i = 0; i < row(); ++i) {
            for(int j = 0; j < getLength(); ++j)
                arr[i][j] = vectors[j][i];
        }
        this->~Matrix();
        length = row();
        vectors = arr;
    }
////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    Matrix reciprocal(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = reciprocal(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix sqrt(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = sqrt(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix factorial(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = factorial(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix ln(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = ln(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix log(const Matrix& m, const Numerical& a) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = log(m[i], a);
        return Matrix(arr, m.getLength());
    }

    Matrix exp(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = exp(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix pow(const Matrix& m, const Numerical& a) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = pow(m[i], a);
        return Matrix(arr, m.getLength());
    }

    Matrix cos(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = cos(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix sin(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = sin(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix tan(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = tan(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix sec(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = sec(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix csc(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = csc(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix cot(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = cot(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arccos(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = arccos(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arcsin(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = arcsin(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arctan(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = arctan(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arcsec(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = arcsec(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arccsc(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = arccsc(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arccot(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = arccot(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix cosh(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = cosh(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix sinh(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = sinh(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix tanh(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = tanh(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix sech(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = sech(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix csch(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = csch(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix coth(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = coth(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arccosh(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = arccosh(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arcsinh(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = arcsinh(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arctanh(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = arctanh(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arcsech(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = arcsech(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arccsch(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = arccsch(m[i]);
        return Matrix(arr, m.getLength());
    }

    Matrix arccoth(const Matrix& m) {
        auto arr = new Vector[m.getLength()];
        for(int i = 0; i < m.getLength(); ++i)
            arr[i] = arccoth(m[i]);
        return Matrix(arr, m.getLength());
    }
}