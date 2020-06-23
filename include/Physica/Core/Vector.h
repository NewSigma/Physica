#ifndef PHYSICA_VECTOR_H
#define PHYSICA_VECTOR_H

#include <iosfwd>
#include "Core/Header/Scalar.h"
#include "Core/Header/Utils/CStyleArray/CStyleArray.h"

namespace Physica::Core {
    class Vector : public CStyleArray<Scalar, Dynamic> {
    public:
        Vector();
        explicit Vector(size_t capacity);
        Vector(size_t length, size_t capacity);
        Vector(const Scalar& n1, const Scalar& n2);
        Vector(const Scalar& n1, const Scalar& n2, const Scalar& n3);
        explicit Vector(const CStyleArray<Scalar, Dynamic>& array);
        explicit Vector(CStyleArray<Scalar, Dynamic>&& array) noexcept;
        Vector(const Vector& vec) = default;
        Vector(Vector&& vec) noexcept;
        ~Vector() = default;
        /* Operators */
        friend std::ostream& operator<<(std::ostream& os, const Vector& v);
        Vector& operator<<(double d) { static_cast<CStyleArray<Scalar, Dynamic>>(*this) << Scalar(d); return *this; }
        Vector& operator=(const Vector& v) noexcept { CStyleArray<Scalar, Dynamic>::operator=(v); return *this; }
        Vector& operator=(Vector&& v) noexcept { CStyleArray<Scalar, Dynamic>::operator=(std::move(v)); return *this; }
        /* Vector Operations */
        [[nodiscard]] Scalar toNorm() const;
        void toUnit();
        /* Getters */
        [[nodiscard]] bool isZeroVector() const;
        [[nodiscard]] Scalar toArg(size_t axe) const;
        /* Helpers */
        static Vector randomVector(size_t length);
        static Vector simplyMultiply(const Vector& v1, const Vector& v2);
    };
    /* Operators */
    Vector operator-(const Vector& v);
    bool operator==(const Vector& v1, const Vector& v2);
    Vector operator+(const Vector& v1, const Vector& v2);
    Vector operator-(const Vector& v1, const Vector& v2);
    Scalar operator*(const Vector& v1, const Vector& v2);
    Vector operator+(const Vector& v, const Scalar& n);
    Vector operator-(const Vector& v, const Scalar& n);
    Vector operator*(const Vector& v, const Scalar& n);
    Vector operator/(const Vector& v, const Scalar& n);
    /* Inline Implements */
    inline void operator+=(Vector& v1, const Vector& v2) { v1 = v1 + v2; }
    inline void operator-=(Vector& v1, const Vector& v2) { v1 = v1 - v2; }
    inline void operator+=(Vector& v, const Scalar& n) { v = v + n; }
    inline void operator-=(Vector& v, const Scalar& n) { v = v - n; }
    inline void operator*=(Vector& v, const Scalar& n) { v = v * n; }
    inline void operator/=(Vector& v, const Scalar& n) { v = v / n; }
////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    Vector reciprocal(const Vector& v);
    Vector sqrt(const Vector& v);
    Vector factorial(const Vector& v);
    Vector ln(const Vector& v);
    Vector log(const Vector& v, const Scalar& a);
    Vector exp(const Vector& v);
    Vector pow(const Vector& v, const Scalar& a);
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
