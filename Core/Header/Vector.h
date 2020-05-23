#ifndef PHYSICA_VECTOR_H
#define PHYSICA_VECTOR_H

#include <iosfwd>
#include "Core/Header/Numerical.h"

namespace Physica::Core {
    class Vector {
        Numerical* numbers;
        size_t length;
    public:
        Vector();
        explicit Vector(size_t length);
        Vector(const Numerical& n1, const Numerical& n2);
        Vector(const Numerical& n1, const Numerical& n2, const Numerical& n3);
        Vector(Numerical*& numbers, size_t length);
        Vector(Numerical*&& numbers, size_t length);
        explicit Vector(const Numerical& arg);
        Vector(const Vector& vec);
        Vector(Vector&& vec) noexcept;
        explicit Vector(const Vector* vector);
        ~Vector();
        //operators
        friend std::ostream& operator<<(std::ostream& os, const Vector& n);
        Numerical& operator[](size_t i) { return numbers[i]; }
        const Numerical& operator[](size_t i) const { return numbers[i]; }
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
        //Vector Operations
        //initVector(int) should be called if and only if numbers equals to nullptr.
        inline void initVector(size_t l) { numbers = new Numerical[l]; length = l; }
        [[nodiscard]] Numerical toNorm() const;
        void toUnit();
        //Public functions
        [[nodiscard]] bool isEmpty() const { return length == 0; }
        [[nodiscard]] size_t getLength() const { return length; };
        [[nodiscard]] bool isZeroVector() const;
        [[nodiscard]] Numerical toArg(size_t axe) const;
    };
    Vector randomVector(size_t length);
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
