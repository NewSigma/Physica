#ifndef PHYSICA_MATRIX_H
#define PHYSICA_MATRIX_H

#include "Vector.h"

namespace Physica::Core {
    class Matrix {
        Vector* vectors;
        size_t length;
    public:
        Matrix();
        Matrix(size_t c, size_t r);
        Matrix(Vector*& vectors, size_t length);
        Matrix(Vector*&& vectors, size_t length);
        Matrix(Matrix& matrix);
        Matrix(Matrix&& matrix) noexcept;
        ~Matrix();
        //Operators
        Vector& operator[](size_t column) { return vectors[column]; }
        const Vector& operator[](size_t column) const { return vectors[column]; }
        const Numerical& operator()(size_t row, size_t column) const { return vectors[column][row]; }
        Matrix& operator=(const Matrix& m) noexcept;
        Matrix& operator=(Matrix&& m) noexcept;
        //Matrix operations
        void rowSwap(size_t r1, size_t r2) noexcept;
        void columnSwap(size_t c1, size_t c2) noexcept;
        void rowEliminate(size_t r1, size_t r2, size_t element);
        void columnEliminate(size_t c1, size_t c2, size_t element);
        void transpose();
        //Public functions
        [[nodiscard]] inline size_t row() const { return vectors[0].getLength(); }
        [[nodiscard]] inline size_t column() const { return length; }
        [[nodiscard]] inline size_t getLength() const { return length; }
        [[nodiscard]] inline bool isEmpty() const { return length == 0; };
    };
    Matrix operator+(const Matrix& m1, const Matrix& m2);
    Matrix operator-(const Matrix& m1, const Matrix& m2);
    Matrix operator*(const Matrix& m1, const Matrix& m2);
    Matrix operator*(const Matrix& m, const Numerical& n);
    void operator+=(Matrix& m1, const Matrix& m2) { m1 = m1 + m2; }
    void operator-=(Matrix& m1, const Matrix& m2) { m1 = m1 - m2; }
    void operator*=(Matrix& m1, const Matrix& m2) { m1 = m1 * m2; }
    void operator*=(Matrix& m, const Numerical& n) { m = m * n; }
////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    Matrix reciprocal(const Matrix& n);
    Matrix sqrt(const Matrix& n);
    Matrix factorial(const Matrix& n);
    Matrix ln(const Matrix& n);
    Matrix log(const Matrix& n, const Numerical& a);
    Matrix exp(const Matrix& n);
    Matrix pow(const Matrix& n, const Numerical& a);
    Matrix cos(const Matrix& n);
    Matrix sin(const Matrix& n);
    Matrix tan(const Matrix& n);
    Matrix sec(const Matrix& n);
    Matrix csc(const Matrix& n);
    Matrix cot(const Matrix& n);
    Matrix arccos(const Matrix& n);
    Matrix arcsin(const Matrix& n);
    Matrix arctan(const Matrix& n);
    Matrix arcsec(const Matrix& n);
    Matrix arccsc(const Matrix& n);
    Matrix arccot(const Matrix& n);
    Matrix cosh(const Matrix& n);
    Matrix sinh(const Matrix& n);
    Matrix tanh(const Matrix& n);
    Matrix sech(const Matrix& n);
    Matrix csch(const Matrix& n);
    Matrix coth(const Matrix& n);
    Matrix arccosh(const Matrix& n);
    Matrix arcsinh(const Matrix& n);
    Matrix arctanh(const Matrix& n);
    Matrix arcsech(const Matrix& n);
    Matrix arccsch(const Matrix& n);
    Matrix arccoth(const Matrix& n);
}

#endif
