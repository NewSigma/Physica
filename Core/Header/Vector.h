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
    int getLength() const { return length; };
    bool isEmpty() const { return length == 0; };
};
////////////////////////////////////////Elementary Functions////////////////////////////////////////////
Vector* reciprocal(const Vector& n);
Vector* sqrt(const Vector& n);
Vector* factorial(const Vector& n);
Vector* ln(const Vector& n);
Vector* log(const Vector& n, const AbstractNum& a);
Vector* exp(const Vector& n);
Vector* pow(const Vector& n, const AbstractNum& a);
Vector* cos(const Vector& n);
Vector* sin(const Vector& n);
Vector* tan(const Vector& n);
Vector* sec(const Vector& n);
Vector* csc(const Vector& n);
Vector* cot(const Vector& n);
Vector* arccos(const Vector& n);
Vector* arcsin(const Vector& n);
Vector* arctan(const Vector& n);
Vector* arcsec(const Vector& n);
Vector* arccsc(const Vector& n);
Vector* arccot(const Vector& n);
Vector* cosh(const Vector& n);
Vector* sinh(const Vector& n);
Vector* tanh(const Vector& n);
Vector* sech(const Vector& n);
Vector* csch(const Vector& n);
Vector* coth(const Vector& n);
Vector* arccosh(const Vector& n);
Vector* arcsinh(const Vector& n);
Vector* arctanh(const Vector& n);
Vector* arcsech(const Vector& n);
Vector* arccsch(const Vector& n);
Vector* arccoth(const Vector& n);

#endif
