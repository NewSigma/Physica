#ifndef PHYSICA_VECTOR_H
#define PHYSICA_VECTOR_H

#include "AbstractNum.h"

class Vector {
public:
    AbstractNum** numbers;
protected:
    int length;
public:
    Vector();
    Vector(AbstractNum* n1, AbstractNum* n2);
    Vector(AbstractNum* n1, AbstractNum* n2, AbstractNum* n3);
    Vector(AbstractNum** numbers, int length);
    Vector(AbstractNum* arg);
    Vector(Vector& vector);
    Vector(Vector* vector);
    ~Vector();

    AbstractNum* toNorm() const;
    AbstractNum* toArg(int index) const;
    void toUnit();

    friend std::ostream& operator<<(std::ostream& os, const Vector& n);
    AbstractNum* operator[](int i) const;
    void operator<<(Vector& v);
    Vector* operator+(Vector& v);
    Vector* operator-(Vector& v);
    AbstractNum* operator*(Vector& v);
    Vector* operator/(Vector& v);
    Vector* operator*(const AbstractNum& n);
    void operator+=(Vector& v);
    void operator-=(Vector& v);
    void operator/=(Vector& v);
    void operator*=(AbstractNum& n);
    Vector* operator-() const;
    bool operator==(Vector& v);
    int getLength() { return length; };
    bool isEmpty() { return length == 0; };
};

#endif
