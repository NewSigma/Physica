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
#ifndef PHYSICA_COMPLEXSCALAR_H
#define PHYSICA_COMPLEXSCALAR_H

#include "Physica/Core/MultiPrecition/Scalar.h"

namespace Physica::Core {
    template<ScalarType type = MultiPrecision, bool errorTrack = true>
    class ComplexScalar {
        Scalar<type, errorTrack> real;
        Scalar<type, errorTrack> imag;
    public:
        ComplexScalar() = default;
        ComplexScalar(Scalar<type, errorTrack> s1, Scalar<type, errorTrack> s2);
        ComplexScalar(const ComplexScalar& c) = default;
        ComplexScalar(ComplexScalar&& c) noexcept;
        /* Operators */
        ComplexScalar& operator=(const ComplexScalar& c);
        ComplexScalar& operator=(ComplexScalar&& c) noexcept;
        void operator<<=(int i) { real <<= i; imag<<= i; }
        void operator>>=(int i) { real >>= i; imag>>= i; }
        bool operator==(const ComplexScalar<type, errorTrack>& c);
        bool operator!=(const ComplexScalar<type, errorTrack>& c) { return !(operator==(c)); }
        bool operator>(const ComplexScalar<type, errorTrack>& c);
        bool operator>=(const ComplexScalar<type, errorTrack>& c) { return !(operator<(c)); }
        bool operator<(const ComplexScalar<type, errorTrack>& c);
        bool operator<=(const ComplexScalar<type, errorTrack>& c) { return !(operator>(c)); }
        /* Helpers */
        void swap(ComplexScalar& c) noexcept;
        static inline ComplexScalar getZero();
        static inline ComplexScalar getOne();
        static inline ComplexScalar getTwo();
        static inline ComplexScalar getRandom();
        /* Getters */
        [[nodiscard]] const Scalar<type, errorTrack>& getReal() const { return real; }
        [[nodiscard]] const Scalar<type, errorTrack>& getImag() const { return imag; }
        [[nodiscard]] bool isZero() { return real.isZero() && imag.isZero(); }
        /* Setters */
        void setReal(const Scalar<type, errorTrack>& s) { real = s; }
        void setImag(const Scalar<type, errorTrack>& s) { imag = s; }
        void setReal(Scalar<type, errorTrack>&& s) { real = std::move(s); }
        void setImag(Scalar<type, errorTrack>&& s) { imag = std::move(s); }
    };
    template<ScalarType type, bool errorTrack>
    [[nodiscard]] inline Scalar<type, errorTrack> norm(const ComplexScalar<type, errorTrack>& c);

    template<ScalarType type, bool errorTrack>
    [[nodiscard]] Scalar<type, errorTrack> arg(const ComplexScalar<type, errorTrack>& c);

    template<ScalarType type, bool errorTrack>
    std::ostream& operator<<(std::ostream& os, const ComplexScalar<type, errorTrack>& c);

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> operator+(
            const ComplexScalar<type, errorTrack>& c1, const ComplexScalar<type, errorTrack>& c2);

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> operator-(
            const ComplexScalar<type, errorTrack>& c1, const ComplexScalar<type, errorTrack>& c2);

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> operator*(
            const ComplexScalar<type, errorTrack>& c1, const ComplexScalar<type, errorTrack>& c2);

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> operator/(
            const ComplexScalar<type, errorTrack>& c1, const ComplexScalar<type, errorTrack>& c2);

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator+(
            const ComplexScalar<type, errorTrack1>& c, const Scalar<type, errorTrack2>& s);

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator-(
            const ComplexScalar<type, errorTrack1>& c, const Scalar<type, errorTrack2>& s);

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator*(
            const ComplexScalar<type, errorTrack1>& c, const Scalar<type, errorTrack2>& s);

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator/(
            const ComplexScalar<type, errorTrack1>& c, const Scalar<type, errorTrack2>& s);

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator+(
            const Scalar<type, errorTrack1>& c, const ComplexScalar<type, errorTrack2>& s);

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator-(
            const Scalar<type, errorTrack1>& c, const ComplexScalar<type, errorTrack2>& s);

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator*(
            const Scalar<type, errorTrack1>& c, const ComplexScalar<type, errorTrack2>& s);

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator/(
            const Scalar<type, errorTrack1>& c, const ComplexScalar<type, errorTrack2>& s);

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> operator<<(const ComplexScalar<type, errorTrack>& c, int i) {
        return ComplexScalar<type, errorTrack>(c.getReal() << i, c.getImag() << i);
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> operator>>(const ComplexScalar<type, errorTrack>& c, int i) {
        return ComplexScalar<type, errorTrack>(c.getReal() >> i, c.getImag() >> i);
    }

    template<ScalarType type, bool errorTrack>
    void swap(ComplexScalar<type, errorTrack>& c1, ComplexScalar<type, errorTrack>& c2) noexcept { c1.swap(c2); }

    template<ScalarType type, bool errorTrack, class T>
    void operator+=(ComplexScalar<type, errorTrack>& c, const T& t) { c = c + t; }

    template<ScalarType type, bool errorTrack, class T>
    void operator-=(ComplexScalar<type, errorTrack>& c, const T& t) { c = c - t; }

    template<ScalarType type, bool errorTrack, class T>
    void operator*=(ComplexScalar<type, errorTrack>& c, const T& t) { c = c * t; }

    template<ScalarType type, bool errorTrack, class T>
    void operator/=(ComplexScalar<type, errorTrack>& c, const T& t) { c = c / t; }
}

#include "ComplexScalarImpl/ComplexScalarImpl.h"

#endif