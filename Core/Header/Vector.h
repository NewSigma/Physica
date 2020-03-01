#ifndef PHYSICA_VECTOR_H
#define PHYSICA_VECTOR_H

#include "RealNumber.h"

class Vector {
public:
    RealNumber** numbers;
private:
    int length;
public:
    Vector();
    Vector(RealNumber** numbers, int length);
    Vector(Vector& vector);
    ~Vector();

    RealNumber* operator[](int i);
    void operator<<(Vector& v);
    Vector* operator+(Vector& v);
    Vector* operator-(Vector& v);
    RealNumber* operator*(Vector& v);
    Vector* operator/(Vector& v);
    Vector* operator*(RealNumber& n);
    void operator+=(Vector& v);
    void operator-=(Vector& v);
    void operator/=(Vector& v);
    void operator*=(RealNumber& n);
    int getLength() { return length; };
    bool isEmpty() { return length == 0; };
};

#endif
