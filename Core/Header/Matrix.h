#ifndef PHYSICA_MATRIX_H
#define PHYSICA_MATRIX_H

#include "Vector.h"

class Matrix {
    Vector* vectors;
    int length;
public:
    Matrix();
    Matrix(Vector* vectors, int length);
    Matrix(Matrix& matrix);
    Matrix(Matrix&& matrix) noexcept;
    ~Matrix();

    Vector& operator[](int index) const;
    Numerical& operator()(int c, int r) const { return vectors[c][r]; }
    Matrix& operator=(const Matrix& m) noexcept;
    Matrix& operator=(Matrix&& m) noexcept;

    inline int row() const { return vectors[0].getLength(); }
    inline int col() const { return length; }
    inline int getLength() const { return length; }
    inline bool isEmpty() const { return length == 0; };
    Matrix inverse() const;
    void transpose();
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

#endif
