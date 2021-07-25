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

#include "Physica/Core/MultiPrecision/Scalar.h"

namespace Physica::Core {
    template<ScalarOption option = MultiPrecision, bool errorTrack = true>
    class ComplexScalar {
        Scalar<option, errorTrack> real;
        Scalar<option, errorTrack> imag;
    public:
        ComplexScalar() = default;
        ComplexScalar(Scalar<option, errorTrack> s1, Scalar<option, errorTrack> s2);
        ComplexScalar(const ComplexScalar& c) = default;
        ComplexScalar(ComplexScalar&& c) noexcept;
        /* Operators */
        ComplexScalar& operator=(const ComplexScalar& c);
        ComplexScalar& operator=(ComplexScalar&& c) noexcept;
        void operator<<=(int i) { real <<= i; imag<<= i; }
        void operator>>=(int i) { real >>= i; imag>>= i; }
        bool operator==(const ComplexScalar<option, errorTrack>& c);
        bool operator!=(const ComplexScalar<option, errorTrack>& c) { return !(operator==(c)); }
        bool operator>(const ComplexScalar<option, errorTrack>& c);
        bool operator>=(const ComplexScalar<option, errorTrack>& c) { return !(operator<(c)); }
        bool operator<(const ComplexScalar<option, errorTrack>& c);
        bool operator<=(const ComplexScalar<option, errorTrack>& c) { return !(operator>(c)); }
        /* Helpers */
        void swap(ComplexScalar& c) noexcept;
        static inline ComplexScalar getZero();
        static inline ComplexScalar getOne();
        static inline ComplexScalar getTwo();
        static inline ComplexScalar getRandom();
        /* Getters */
        [[nodiscard]] const Scalar<option, errorTrack>& getReal() const { return real; }
        [[nodiscard]] const Scalar<option, errorTrack>& getImag() const { return imag; }
        [[nodiscard]] bool isZero() { return real.isZero() && imag.isZero(); }
        /* Setters */
        void setReal(const Scalar<option, errorTrack>& s) { real = s; }
        void setImag(const Scalar<option, errorTrack>& s) { imag = s; }
        void setReal(Scalar<option, errorTrack>&& s) { real = std::move(s); }
        void setImag(Scalar<option, errorTrack>&& s) { imag = std::move(s); }
    };
    template<ScalarOption option, bool errorTrack>
    [[nodiscard]] inline Scalar<option, errorTrack> norm(const ComplexScalar<option, errorTrack>& c);

    template<ScalarOption option, bool errorTrack>
    [[nodiscard]] Scalar<option, errorTrack> arg(const ComplexScalar<option, errorTrack>& c);

    template<ScalarOption option, bool errorTrack>
    std::ostream& operator<<(std::ostream& os, const ComplexScalar<option, errorTrack>& c);

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> operator+(
            const ComplexScalar<option, errorTrack>& c1, const ComplexScalar<option, errorTrack>& c2);

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> operator-(
            const ComplexScalar<option, errorTrack>& c1, const ComplexScalar<option, errorTrack>& c2);

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> operator*(
            const ComplexScalar<option, errorTrack>& c1, const ComplexScalar<option, errorTrack>& c2);

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> operator/(
            const ComplexScalar<option, errorTrack>& c1, const ComplexScalar<option, errorTrack>& c2);

    template<ScalarOption option, bool errorTrack1, bool errorTrack2>
    ComplexScalar<option, errorTrack1 || errorTrack2> operator+(
            const ComplexScalar<option, errorTrack1>& c, const Scalar<option, errorTrack2>& s);

    template<ScalarOption option, bool errorTrack1, bool errorTrack2>
    ComplexScalar<option, errorTrack1 || errorTrack2> operator-(
            const ComplexScalar<option, errorTrack1>& c, const Scalar<option, errorTrack2>& s);

    template<ScalarOption option, bool errorTrack1, bool errorTrack2>
    ComplexScalar<option, errorTrack1 || errorTrack2> operator*(
            const ComplexScalar<option, errorTrack1>& c, const Scalar<option, errorTrack2>& s);

    template<ScalarOption option, bool errorTrack1, bool errorTrack2>
    ComplexScalar<option, errorTrack1 || errorTrack2> operator/(
            const ComplexScalar<option, errorTrack1>& c, const Scalar<option, errorTrack2>& s);

    template<ScalarOption option, bool errorTrack1, bool errorTrack2>
    ComplexScalar<option, errorTrack1 || errorTrack2> operator+(
            const Scalar<option, errorTrack1>& c, const ComplexScalar<option, errorTrack2>& s);

    template<ScalarOption option, bool errorTrack1, bool errorTrack2>
    ComplexScalar<option, errorTrack1 || errorTrack2> operator-(
            const Scalar<option, errorTrack1>& c, const ComplexScalar<option, errorTrack2>& s);

    template<ScalarOption option, bool errorTrack1, bool errorTrack2>
    ComplexScalar<option, errorTrack1 || errorTrack2> operator*(
            const Scalar<option, errorTrack1>& c, const ComplexScalar<option, errorTrack2>& s);

    template<ScalarOption option, bool errorTrack1, bool errorTrack2>
    ComplexScalar<option, errorTrack1 || errorTrack2> operator/(
            const Scalar<option, errorTrack1>& c, const ComplexScalar<option, errorTrack2>& s);

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> operator<<(const ComplexScalar<option, errorTrack>& c, int i) {
        return ComplexScalar<option, errorTrack>(c.getReal() << i, c.getImag() << i);
    }

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> operator>>(const ComplexScalar<option, errorTrack>& c, int i) {
        return ComplexScalar<option, errorTrack>(c.getReal() >> i, c.getImag() >> i);
    }

    template<ScalarOption option, bool errorTrack>
    void swap(ComplexScalar<option, errorTrack>& c1, ComplexScalar<option, errorTrack>& c2) noexcept { c1.swap(c2); }

    template<ScalarOption option, bool errorTrack, class T>
    void operator+=(ComplexScalar<option, errorTrack>& c, const T& t) { c = c + t; }

    template<ScalarOption option, bool errorTrack, class T>
    void operator-=(ComplexScalar<option, errorTrack>& c, const T& t) { c = c - t; }

    template<ScalarOption option, bool errorTrack, class T>
    void operator*=(ComplexScalar<option, errorTrack>& c, const T& t) { c = c * t; }

    template<ScalarOption option, bool errorTrack, class T>
    void operator/=(ComplexScalar<option, errorTrack>& c, const T& t) { c = c / t; }
}

#include "ComplexScalarImpl/ComplexScalarImpl.h"

#endif