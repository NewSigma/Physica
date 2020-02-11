#ifndef PHYSICA_VECTOR_H
#define PHYSICA_VECTOR_H

#include "RealNumber.h"
#include "Matrix.h"

class Vector {
private:
    RealNumber** numbers;
    int length;
public:
    Vector();
    Vector(RealNumber** numbers, int length);
    ~Vector();

    RealNumber* operator[](int i);
    Vector* operator+(Vector& v);
    Vector* operator-(Vector& v);
    RealNumber* operator*(Vector& v);
    Vector* operator/(Vector& v);
    int getLength() { return length; };
    bool isEmpty() { return length == 0; };
//////////////////////////////friends///////////////////////////////////
    friend void Matrix::toColMatrix();
    friend void Matrix::toRowMatrix();
};

#endif
