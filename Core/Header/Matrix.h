#ifndef PHYSICA_MATRIX_H
#define PHYSICA_MATRIX_H

#include "AbstractNum.h"

class Vector;

class Matrix {
    //When Type = COLUMN, vectors stores column vectors.
    enum Type {
        ROW,
        COLUMN
    };
public:
    Matrix();
    Matrix(Vector** vectors, int length, Type type = COLUMN);
    Matrix(Matrix& matrix);
    ~Matrix();

    Vector* operator[](int index) const;
    void operator<<(Matrix& m);
    Matrix* operator+(Matrix& m);
    Matrix* operator-(Matrix& m);
    Matrix* operator*(Matrix& m);
    Matrix* operator/(Matrix& m);
    Matrix* operator*(AbstractNum& n);
    void operator+=(Matrix& m);
    void operator-=(Matrix& m);
    void operator*=(Matrix& m);
    void operator/=(Matrix& m);
    void operator*=(AbstractNum& n);

    int getRow() const { return row; };
    int getCol() const { return col; };
    int getLength() const { if(type) return row; else return col; };
    Type getType() const { return type; }
    inline bool isEmpty() const { return row == 0; };
    AbstractNum* get(int c, int r);
    void toColMatrix();
    void toRowMatrix();
    void changeType();
    Matrix* inverse();
    void transpose();
private:
    Vector** vectors;
    int col, row;
    //if(type) <==> if(type == COLUMN)
    Type type;
};
////////////////////////////////////////Elementary Functions////////////////////////////////////////////
Matrix* reciprocal(const Matrix& n);
Matrix* sqrt(const Matrix& n);
Matrix* factorial(const Matrix& n);
Matrix* ln(const Matrix& n);
Matrix* log(const Matrix& n, const AbstractNum& a);
Matrix* exp(const Matrix& n);
Matrix* pow(const Matrix& n, const AbstractNum& a);
Matrix* cos(const Matrix& n);
Matrix* sin(const Matrix& n);
Matrix* tan(const Matrix& n);
Matrix* sec(const Matrix& n);
Matrix* csc(const Matrix& n);
Matrix* cot(const Matrix& n);
Matrix* arccos(const Matrix& n);
Matrix* arcsin(const Matrix& n);
Matrix* arctan(const Matrix& n);
Matrix* arcsec(const Matrix& n);
Matrix* arccsc(const Matrix& n);
Matrix* arccot(const Matrix& n);
Matrix* cosh(const Matrix& n);
Matrix* sinh(const Matrix& n);
Matrix* tanh(const Matrix& n);
Matrix* sech(const Matrix& n);
Matrix* csch(const Matrix& n);
Matrix* coth(const Matrix& n);
Matrix* arccosh(const Matrix& n);
Matrix* arcsinh(const Matrix& n);
Matrix* arctanh(const Matrix& n);
Matrix* arcsech(const Matrix& n);
Matrix* arccsch(const Matrix& n);
Matrix* arccoth(const Matrix& n);

#endif
