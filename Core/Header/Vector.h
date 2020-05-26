#ifndef PHYSICA_VECTOR_H
#define PHYSICA_VECTOR_H

#include <iosfwd>
#include "Core/Header/Numerical.h"

namespace Physica::Core {
    class Vector {
        Numerical* __restrict numbers;
        size_t length;
        size_t capacity;
    public:
        Vector();
        explicit Vector(size_t length);
        Vector(const Numerical& n1, const Numerical& n2);
        Vector(const Numerical& n1, const Numerical& n2, const Numerical& n3);
        Vector(Numerical*& numbers, size_t length);
        Vector(Numerical*&& numbers, size_t length);
        Vector(const Vector& vec);
        Vector(Vector&& vec) noexcept;
        ~Vector();
        /* Operators */
        friend std::ostream& operator<<(std::ostream& os, const Vector& n);
        Numerical& operator[](size_t i) { return numbers[i]; }
        const Numerical& operator[](size_t i) const { return numbers[i]; }
        Vector& operator=(const Vector& v) noexcept;
        Vector& operator=(Vector&& v) noexcept;
        Vector operator/(const Vector& v) const;
        /* Vector Operations */
        void initVector(size_t l) { delete[] numbers; numbers = new Numerical[l]; length = l; }
        [[nodiscard]] Numerical toNorm() const;
        void toUnit();
        void resize(size_t size);
        void squeeze();
        void append(Numerical n) noexcept;
        void append(Vector v) noexcept;
        Vector cut(size_t from);
        Numerical cutLast();
        /* Getters */
        void swap(Vector& v) noexcept;
        [[nodiscard]] Vector subVector(size_t from, size_t to);
        /*
         * Another version of subVector(size_t, size_t). From a index to the end of vector.
         */
        [[nodiscard]] Vector subVector(size_t from) { return subVector(from, length); }
        [[nodiscard]] bool isEmpty() const { return length == 0; }
        [[nodiscard]] size_t getLength() const { return length; };
        [[nodiscard]] size_t getCapacity() const { return capacity; }
        [[nodiscard]] bool isZeroVector() const;
        [[nodiscard]] Numerical toArg(size_t axe) const;
    };
    Vector randomVector(size_t length);
    Vector simplyMultiply(const Vector& v1, const Vector& v2);
    /* Operators */
    Vector operator-(const Vector& v);
    bool operator==(const Vector& v1, const Vector& v2);
    Vector operator+(const Vector& v1, const Vector& v2);
    Vector operator-(const Vector& v1, const Vector& v2);
    Numerical operator*(const Vector& v1, const Vector& v2);
    Vector operator+(const Vector& v, const Numerical& n);
    Vector operator-(const Vector& v, const Numerical& n);
    Vector operator*(const Vector& v, const Numerical& n);
    Vector operator/(const Vector& v, const Numerical& n);
    /* Inline Implements */
    inline void swap(Vector& v1, Vector& v2) noexcept { v1.swap(v2); }
    inline void operator+=(Vector& v1, const Vector& v2) { v1 = v1 + v2; };
    inline void operator-=(Vector& v1, const Vector& v2) { v1 = v1 - v2; };
    inline void operator+=(Vector& v, const Numerical& n) { v = v + n; };
    inline void operator-=(Vector& v, const Numerical& n) { v = v - n; };
    inline void operator*=(Vector& v, const Numerical& n) { v = v * n; };
    inline void operator/=(Vector& v, const Numerical& n) { v = v / n; };
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
