#ifndef PHYSICA_VECTOR_H
#define PHYSICA_VECTOR_H

#include "RealNumber.h"

class Vector {
public:
    RealNumber** numbers;
    int length;

    Vector(RealNumber** numbers, int length);
    ~Vector();

    Vector& operator+(Vector& v);
    Vector& operator-(Vector& v);
    Vector& operator*(Vector& v);
    Vector& operator/(Vector& v);
};

#endif
