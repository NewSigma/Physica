#ifndef PHYSICA_MATRIX_H
#define PHYSICA_MATRIX_H

#include "AbstractNum.h"

class Vector;

class Matrix {
private:
    //When Type = COLUMN, vectors stores column vectors.
    enum Type {
        ROW,
        COLUMN
    };

    Vector** vectors;
    int col, row;
    //if(type) <==> if(type == COLUMN)
    Type type;
public:
    Matrix();
    Matrix(Vector** vectors, int length, Type type = COLUMN);
    Matrix(Matrix& matrix);
    ~Matrix();

    Vector* operator[](int index);
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

    int getRow() { return row; };
    int getCol() { return col; };
    int getLength() { if(type) return row; else return col; };
    bool isEmpty() { return row == 0; };
    AbstractNum* get(int c, int r);
    void toColMatrix();
    void toRowMatrix();
    void changeType();
    Matrix* inverse();
    void transpose();
};

#endif
