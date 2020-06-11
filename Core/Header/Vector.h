#ifndef PHYSICA_VECTOR_H
#define PHYSICA_VECTOR_H

#include <iosfwd>
#include "Core/Header/Numerical.h"

namespace Physica::Core {
    /*
     * TODO:
     * 1.new class MatrixPrivate, VectorPrivate (Handle resource, template, malloc, placement new)
     * 2.Constructor Matrix(MatrixPrivate) Vector(VectorPrivate)
     * 3.remove grow()
     * 4.add setLength()
     */
    class Vector {
        Numerical* __restrict numbers;
        size_t length;
        size_t capacity;
    public:
        Vector();
        explicit Vector(size_t length);
        Vector(const Numerical& n1, const Numerical& n2);
        Vector(const Numerical& n1, const Numerical& n2, const Numerical& n3);
        Vector(const Vector& vec);
        Vector(Vector&& vec) noexcept;
        ~Vector();
        /* Operators */
        friend std::ostream& operator<<(std::ostream& os, const Vector& v);
        Vector& operator<<(Numerical n);
        Numerical& operator[](size_t i) { Q_ASSERT(i < length); return numbers[i]; }
        const Numerical& operator[](size_t i) const { Q_ASSERT(i < length); return numbers[i]; }
        Vector& operator=(const Vector& v) noexcept;
        Vector& operator=(Vector&& v) noexcept;
        Vector operator/(const Vector& v) const;
        /* Vector Operations */
        [[nodiscard]] Numerical toNorm() const;
        void toUnit();
        /* Getters */
        [[nodiscard]] Vector subVector(size_t from, size_t to);
        // Another version of subVector(size_t, size_t). From a index to the end of vector.
        [[nodiscard]] Vector subVector(size_t from) { return subVector(from, length); }
        [[nodiscard]] bool isEmpty() const { return length == 0; }
        [[nodiscard]] size_t getLength() const { return length; };
        [[nodiscard]] size_t getCapacity() const { return capacity; }
        [[nodiscard]] bool isZeroVector() const;
        [[nodiscard]] Numerical toArg(size_t axe) const;
        /* Helpers */
        void resize(size_t size);
        void squeeze();
        void append(Numerical n) noexcept;
        void append(Vector v) noexcept;
        Vector cut(size_t from);
        Numerical cutLast();
        inline void grow(Numerical n);
        void swap(Vector& v) noexcept;
        static Vector randomVector(size_t length);
        static Vector simplyMultiply(const Vector& v1, const Vector& v2);
    protected:
        Vector(Numerical*& numbers, size_t length);
        Vector(Numerical*&& numbers, size_t length);
    };
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
    //A low configuration version of append(Numerical n), which do not adjust the allocated memory.
    inline void Vector::grow(Numerical n) {
        Q_ASSERT(length < capacity);
        new (numbers + length) Numerical(std::move(n));
        ++length;
    }
    inline void swap(Vector& v1, Vector& v2) noexcept { v1.swap(v2); }
    inline void operator+=(Vector& v1, const Vector& v2) { v1 = v1 + v2; };
    inline void operator-=(Vector& v1, const Vector& v2) { v1 = v1 - v2; };
    inline void operator+=(Vector& v, const Numerical& n) { v = v + n; };
    inline void operator-=(Vector& v, const Numerical& n) { v = v - n; };
    inline void operator*=(Vector& v, const Numerical& n) { v = v * n; };
    inline void operator/=(Vector& v, const Numerical& n) { v = v / n; };
////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    Vector reciprocal(const Vector& v);
    Vector sqrt(const Vector& v);
    Vector factorial(const Vector& v);
    Vector ln(const Vector& v);
    Vector log(const Vector& v, const Numerical& a);
    Vector exp(const Vector& v);
    Vector pow(const Vector& v, const Numerical& a);
    Vector cos(const Vector& v);
    Vector sin(const Vector& v);
    Vector tan(const Vector& v);
    Vector sec(const Vector& v);
    Vector csc(const Vector& v);
    Vector cot(const Vector& v);
    Vector arccos(const Vector& v);
    Vector arcsin(const Vector& v);
    Vector arctan(const Vector& v);
    Vector arcsec(const Vector& v);
    Vector arccsc(const Vector& v);
    Vector arccot(const Vector& v);
    Vector cosh(const Vector& v);
    Vector sinh(const Vector& v);
    Vector tanh(const Vector& v);
    Vector sech(const Vector& v);
    Vector csch(const Vector& v);
    Vector coth(const Vector& v);
    Vector arccosh(const Vector& v);
    Vector arcsinh(const Vector& v);
    Vector arctanh(const Vector& v);
    Vector arcsech(const Vector& v);
    Vector arccsch(const Vector& v);
    Vector arccoth(const Vector& v);
}

#endif
