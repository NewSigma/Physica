/*
 * Copyright 2020-2021 WeiBo He.
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
#pragma once

#include "Physica/Core/MultiPrecision/Scalar.h"

namespace Physica::Core {
    template<class AnyScalar> class ComplexScalar;

    namespace Internal {
        template<class AnyScalar>
        class Traits<ComplexScalar<AnyScalar>> {
        public:
            static constexpr ScalarOption option = Traits<AnyScalar>::option;
            static constexpr bool errorTrack = Traits<AnyScalar>::errorTrack;
            static constexpr bool isComplex = true;
        };
    }

    template<class AnyScalar>
    class ComplexScalar : public ScalarBase<ComplexScalar<AnyScalar>> {
    public:
        using ScalarType = ComplexScalar<AnyScalar>;
    private:
        AnyScalar real;
        AnyScalar imag;
    public:
        ComplexScalar() = default;
        ComplexScalar(const ScalarBase<AnyScalar>& real_);
        ComplexScalar(const ScalarBase<AnyScalar>& real_, const ScalarBase<AnyScalar>& imag_);
        ComplexScalar(std::initializer_list<AnyScalar> list);
        ComplexScalar(const ComplexScalar& c) = default;
        ComplexScalar(ComplexScalar&& c) noexcept;
        /* Operators */
        ComplexScalar& operator=(const ComplexScalar& c);
        ComplexScalar& operator=(ComplexScalar&& c) noexcept;
        ComplexScalar& operator=(const ScalarBase<AnyScalar>& s);
        void operator<<=(int i) { real <<= i; imag<<= i; }
        void operator>>=(int i) { real >>= i; imag>>= i; }
        bool operator==(const ComplexScalar<AnyScalar>& c);
        bool operator!=(const ComplexScalar<AnyScalar>& c) { return !(operator==(c)); }
        bool operator>(const ComplexScalar<AnyScalar>& c);
        bool operator>=(const ComplexScalar<AnyScalar>& c) { return !(operator<(c)); }
        bool operator<(const ComplexScalar<AnyScalar>& c);
        bool operator<=(const ComplexScalar<AnyScalar>& c) { return !(operator>(c)); }
        /* Helpers */
        void swap(ComplexScalar& c) noexcept;
        static inline ComplexScalar Zero();
        static inline ComplexScalar One();
        static inline ComplexScalar Two();
        static inline ComplexScalar Random();
        /* Getters */
        [[nodiscard]] inline AnyScalar norm();
        [[nodiscard]] AnyScalar arg();
        [[nodiscard]] const AnyScalar& getReal() const { return real; }
        [[nodiscard]] const AnyScalar& getImag() const { return imag; }
        [[nodiscard]] bool isZero() { return real.isZero() && imag.isZero(); }
        /* Setters */
        void setReal(const AnyScalar& s) { real = s; }
        void setImag(const AnyScalar& s) { imag = s; }
    };
    template<class AnyScalar>
    std::ostream& operator<<(std::ostream& os, const ComplexScalar<AnyScalar>& c);

    template<class AnyScalar>
    ComplexScalar<AnyScalar> operator+(
            const ComplexScalar<AnyScalar>& c1, const ComplexScalar<AnyScalar>& c2);

    template<class AnyScalar>
    ComplexScalar<AnyScalar> operator-(
            const ComplexScalar<AnyScalar>& c1, const ComplexScalar<AnyScalar>& c2);

    template<class AnyScalar>
    ComplexScalar<AnyScalar> operator*(
            const ComplexScalar<AnyScalar>& c1, const ComplexScalar<AnyScalar>& c2);

    template<class AnyScalar>
    ComplexScalar<AnyScalar> operator/(
            const ComplexScalar<AnyScalar>& c1, const ComplexScalar<AnyScalar>& c2);

    template<class AnyScalar1, class AnyScalar2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type> operator+(
            const ComplexScalar<AnyScalar1>& c, const ScalarBase<AnyScalar2>& s);

    template<class AnyScalar1, class AnyScalar2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type> operator-(
            const ComplexScalar<AnyScalar1>& c, const ScalarBase<AnyScalar2>& s);

    template<class AnyScalar1, class AnyScalar2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type> operator*(
            const ComplexScalar<AnyScalar1>& c, const ScalarBase<AnyScalar2>& s);

    template<class AnyScalar1, class AnyScalar2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type> operator/(
            const ComplexScalar<AnyScalar1>& c, const ScalarBase<AnyScalar2>& s);

    template<class AnyScalar1, class AnyScalar2>
    inline ComplexScalar<typename Internal::BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type> operator+(
            const ScalarBase<AnyScalar1>& s, const ComplexScalar<AnyScalar2>& c);

    template<class AnyScalar1, class AnyScalar2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type> operator-(
            const ScalarBase<AnyScalar1>& s, const ComplexScalar<AnyScalar2>& c);

    template<class AnyScalar1, class AnyScalar2>
    inline ComplexScalar<typename Internal::BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type> operator*(
            const ScalarBase<AnyScalar1>& s, const ComplexScalar<AnyScalar2>& c);

    template<class AnyScalar1, class AnyScalar2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type> operator/(
            const ScalarBase<AnyScalar1>& s, const ComplexScalar<AnyScalar2>& c);

    template<class AnyScalar>
    ComplexScalar<AnyScalar> operator<<(const ComplexScalar<AnyScalar>& c, int i) {
        return ComplexScalar<AnyScalar>(c.getReal() << i, c.getImag() << i);
    }

    template<class AnyScalar>
    ComplexScalar<AnyScalar> operator>>(const ComplexScalar<AnyScalar>& c, int i) {
        return ComplexScalar<AnyScalar>(c.getReal() >> i, c.getImag() >> i);
    }

    template<class AnyScalar>
    void swap(ComplexScalar<AnyScalar>& c1, ComplexScalar<AnyScalar>& c2) noexcept { c1.swap(c2); }

    template<class AnyScalar, class T>
    void operator+=(ComplexScalar<AnyScalar>& c, const T& t) { c = c + t; }

    template<class AnyScalar, class T>
    void operator-=(ComplexScalar<AnyScalar>& c, const T& t) { c = c - t; }

    template<class AnyScalar, class T>
    void operator*=(ComplexScalar<AnyScalar>& c, const T& t) { c = c * t; }

    template<class AnyScalar, class T>
    void operator/=(ComplexScalar<AnyScalar>& c, const T& t) { c = c / t; }
}

#include "ComplexScalarImpl/ComplexScalarImpl.h"
