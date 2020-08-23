/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef PHYSICA_VECTOR_H
#define PHYSICA_VECTOR_H

#include <iosfwd>
#include "Physica/Core/MultiPrecition/Scalar.h"
#include "Physica/Core/Utils/Container/CStyleArray/CStyleArray.h"

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
        Vector(std::initializer_list<T> list);
        Vector(const Vector<T, maxLength>& vec);
        Vector(Vector<T, maxLength>&& vec) noexcept;
        ~Vector() = default;
        /* Operators */
        Vector<T, maxLength>& operator=(const Vector<T, maxLength>& v) noexcept { Base::operator=(v); return *this; }
        Vector<T, maxLength>& operator=(Vector<T, maxLength>&& v) noexcept { Base::operator=(std::move(v)); return *this; }
        /* Vector Operations */
        Vector& toOpposite();
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
