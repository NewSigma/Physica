#ifndef PHYSICA_VECTOR_H
#define PHYSICA_VECTOR_H

#include <iosfwd>
#include "Physica/Core/MultiPrecition/Scalar.h"
#include "Physica/Core/Utils/CStyleArray/CStyleArray.h"

namespace Physica::Core {
    /*!
     * T must be either Scalar or ComplexScalar.
     */
    template<class T = MultiScalar, size_t maxLength = Dynamic>
    class Vector : public CStyleArray<T, maxLength> {
        typedef CStyleArray<T, maxLength> Base;
    public:
        Vector();
        explicit Vector(size_t length);
        explicit Vector(const Base& array);
        explicit Vector(Base&& array) noexcept;
        Vector(const Vector<T, maxLength>& vec);
        Vector(Vector<T, maxLength>&& vec) noexcept;
        ~Vector() = default;
        /* Operators */
        Vector<T, maxLength>& operator=(const Vector<T, maxLength>& v) noexcept { Base::operator=(v); return *this; }
        Vector<T, maxLength>& operator=(Vector<T, maxLength>&& v) noexcept { Base::operator=(std::move(v)); return *this; }
        /* Vector Operations */
        [[nodiscard]] T toNorm() const;
        void toUnit();
        /* Getters */
        [[nodiscard]] bool isZero() const;
        /* Helpers */
        static Vector<T, Dynamic> zeroVector(size_t len);
        static Vector<T, Dynamic> randomVector(size_t len);
        static Vector simplyMultiply(const Vector<T, maxLength>& v1, const Vector<T, maxLength>& v2);
    };
    /* Operators */
    template<class T, size_t maxLength>
    std::ostream& operator<<(std::ostream& os, const Vector<T, maxLength>& v);

    template<class T, size_t maxLength>
    Vector<T, maxLength> operator-(const Vector<T, maxLength>& v);

    template<class T, size_t maxLength>
    Vector<T, maxLength> operator+(const Vector<T, maxLength>& v1, const Vector<T, maxLength>& v2);

    template<class T, size_t maxLength>
    Vector<T, maxLength> operator-(const Vector<T, maxLength>& v1, const Vector<T, maxLength>& v2);

    template<class T, size_t maxLength>
    T operator*(const Vector<T, maxLength>& v1, const Vector<T, maxLength>& v2);

    template<class T, size_t maxLength>
    Vector<T, maxLength> operator+(const Vector<T, maxLength>& v, const T& n);

    template<class T, size_t maxLength>
    Vector<T, maxLength> operator-(const Vector<T, maxLength>& v, const T& n);

    template<class T, size_t maxLength>
    Vector<T, maxLength> operator*(const Vector<T, maxLength>& v, const T& n);

    template<class T, size_t maxLength>
    Vector<T, maxLength> operator/(const Vector<T, maxLength>& v, const T& n);
    /* Inline Implements */
    template<class T, size_t maxLength>
    inline void operator+=(Vector<T, maxLength>& v1, const Vector<T, maxLength>& v2) { v1 = v1 + v2; }

    template<class T, size_t maxLength>
    inline void operator-=(Vector<T, maxLength>& v1, const Vector<T, maxLength>& v2) { v1 = v1 - v2; }

    template<class T, size_t maxLength>
    inline void operator+=(Vector<T, maxLength>& v, const T& n) { v = v + n; }

    template<class T, size_t maxLength>
    inline void operator-=(Vector<T, maxLength>& v, const T& n) { v = v - n; }

    template<class T, size_t maxLength>
    inline void operator*=(Vector<T, maxLength>& v, const T& n) { v = v * n; }

    template<class T, size_t maxLength>
    inline void operator/=(Vector<T, maxLength>& v, const T& n) { v = v / n; }
}

#include "VectorImpl.h"

#endif
