#ifndef PHYSICA_VECTOR_H
#define PHYSICA_VECTOR_H

#include <iosfwd>
#include "Core/Header/Numerical.h"

namespace Physica::Core {
    class Vector {
        Numerical* numbers;
        int length;
    public:
        Vector();
        explicit Vector(int length);
        Vector(const Numerical& n1, const Numerical& n2);
        Vector(const Numerical& n1, const Numerical& n2, const Numerical& n3);
        Vector(Numerical*& numbers, int length);
        Vector(Numerical*&& numbers, int length);
        explicit Vector(const Numerical& arg);
        Vector(const Vector& vec);
        Vector(Vector&& vec) noexcept;
        explicit Vector(const Vector* vector);
        ~Vector();

        Numerical toNorm() const;
        void toUnit();

        friend std::ostream& operator<<(std::ostream& os, const Vector& n);
        Numerical& operator[](int i) { return numbers[i]; }
        const Numerical& operator[](int i) const { return numbers[i]; }
        Vector& operator=(const Vector& v) noexcept;
        Vector& operator=(Vector&& v) noexcept;
        Vector operator+(const Vector& v) const;
        Vector operator-(const Vector& v) const;
        Vector operator*(const Numerical& n) const;
        Numerical operator*(const Vector& v) const;
        Vector operator/(const Vector& v) const;
        void operator+=(const Vector& v) { *this = *this + v; };
        void operator-=(const Vector& v) { *this = *this - v; };
        void operator*=(Numerical& n) { *this = *this * n; };
        Vector operator-() const;
        bool operator==(const Vector& v) const;

        bool isEmpty() const { return length == 0; }
        int getLength() const { return length; };
        bool isZeroVector() const;
        Numerical toArg(int axe) const;
    };
    Vector randomVector(int length);
////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    Vector reciprocal(const Vector& n);
    Vector sqrt(const Vector& n);
    Vector factorial(const Vector& n);
    Vector ln(const Vector& n);
    Vector log(const Vector& n, const Numerical& a);
    Vector exp(const Vector& n);
    Vector pow(const Vector& n, const Numerical& a);
    Vector cos(const Vector& n);
    Vector sin(const Vector& n);
    Vector tan(const Vector& n);
    Vector sec(const Vector& n);
    Vector csc(const Vector& n);
    Vector cot(const Vector& n);
    Vector arccos(const Vector& n);
    Vector arcsin(const Vector& n);
    Vector arctan(const Vector& n);
    Vector arcsec(const Vector& n);
    Vector arccsc(const Vector& n);
    Vector arccot(const Vector& n);
    Vector cosh(const Vector& n);
    Vector sinh(const Vector& n);
    Vector tanh(const Vector& n);
    Vector sech(const Vector& n);
    Vector csch(const Vector& n);
    Vector coth(const Vector& n);
    Vector arccosh(const Vector& n);
    Vector arcsinh(const Vector& n);
    Vector arctanh(const Vector& n);
    Vector arcsech(const Vector& n);
    Vector arccsch(const Vector& n);
    Vector arccoth(const Vector& n);
}

#endif
