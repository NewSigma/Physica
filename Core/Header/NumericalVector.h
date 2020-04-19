#ifndef PHYSICA_NUMERICALVECTOR_H
#define PHYSICA_NUMERICALVECTOR_H

#include <iosfwd>

class Numerical;
class Node;

class NumericalVector {
    Numerical** numbers;
    int length;
public:
    NumericalVector(unsigned int l);
    NumericalVector(const Numerical& n1, const Numerical& n2);
    NumericalVector(const Numerical& n1, const Numerical& n2, const Numerical& n3);
    NumericalVector(Numerical** numbers, int length);
    explicit NumericalVector(const Numerical& arg);
    NumericalVector(const NumericalVector& vector);
    NumericalVector(NumericalVector&& NumericalVector) noexcept;
    explicit NumericalVector(const NumericalVector* vector);
    ~NumericalVector();

    Numerical toNorm() const;
    void toUnit();

    friend std::ostream& operator<<(std::ostream& os, const NumericalVector& n);
    Numerical& operator[](int i) const;
    NumericalVector& operator=(const NumericalVector& v) noexcept;
    NumericalVector& operator=(NumericalVector&& v) noexcept;
    NumericalVector operator+(const NumericalVector& v) const;
    NumericalVector operator-(const NumericalVector& v) const;
    NumericalVector operator*(const Numerical& n) const;
    Numerical operator*(const NumericalVector& v) const;
    NumericalVector operator/(const NumericalVector& v) const;
    void operator+=(const NumericalVector& v) { *this = *this + v; };
    void operator-=(const NumericalVector& v) { *this = *this - v; };
    void operator*=(Numerical& n) { *this = *this * n; };
    NumericalVector operator-() const;
    bool operator==(const NumericalVector& v) const;
    int getLength() const { return length; };
    bool isEmpty() const { return length == 0; };
};
////////////////////////////////////////Elementary Functions////////////////////////////////////////////
NumericalVector reciprocal(const NumericalVector& n);
NumericalVector sqrt(const NumericalVector& n);
NumericalVector factorial(const NumericalVector& n);
NumericalVector ln(const NumericalVector& n);
NumericalVector log(const NumericalVector& n, const Numerical& a);
NumericalVector exp(const NumericalVector& n);
NumericalVector pow(const NumericalVector& n, const Numerical& a);
NumericalVector cos(const NumericalVector& n);
NumericalVector sin(const NumericalVector& n);
NumericalVector tan(const NumericalVector& n);
NumericalVector sec(const NumericalVector& n);
NumericalVector csc(const NumericalVector& n);
NumericalVector cot(const NumericalVector& n);
NumericalVector arccos(const NumericalVector& n);
NumericalVector arcsin(const NumericalVector& n);
NumericalVector arctan(const NumericalVector& n);
NumericalVector arcsec(const NumericalVector& n);
NumericalVector arccsc(const NumericalVector& n);
NumericalVector arccot(const NumericalVector& n);
NumericalVector cosh(const NumericalVector& n);
NumericalVector sinh(const NumericalVector& n);
NumericalVector tanh(const NumericalVector& n);
NumericalVector sech(const NumericalVector& n);
NumericalVector csch(const NumericalVector& n);
NumericalVector coth(const NumericalVector& n);
NumericalVector arccosh(const NumericalVector& n);
NumericalVector arcsinh(const NumericalVector& n);
NumericalVector arctanh(const NumericalVector& n);
NumericalVector arcsech(const NumericalVector& n);
NumericalVector arccsch(const NumericalVector& n);
NumericalVector arccoth(const NumericalVector& n);

#endif
