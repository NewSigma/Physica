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
    template<class ScalarType> class ComplexScalar;

    namespace Internal {
        template<class ScalarType>
        class Traits<ComplexScalar<ScalarType>> {
        public:
            using RealType = ScalarType;
            using ComplexType = ComplexScalar<ScalarType>;
            static constexpr ScalarOption option = Traits<ScalarType>::option;
            static constexpr bool errorTrack = Traits<ScalarType>::errorTrack;
            static constexpr bool isComplex = true;
        };
    }

    template<class ScalarType>
    class ComplexScalar : public ScalarBase<ComplexScalar<ScalarType>> {
        static_assert(!ScalarType::isComplex);
    private:
        ScalarType real;
        ScalarType imag;
    public:
        ComplexScalar() = default;
        ComplexScalar(const ScalarType& real_);
        ComplexScalar(const ScalarType& real_, const ScalarType& imag_);
        ComplexScalar(std::initializer_list<ScalarType> list);
        ComplexScalar(const ComplexScalar& c) = default;
        ComplexScalar(ComplexScalar&& c) noexcept;
        /* Operators */
        ComplexScalar& operator=(const ComplexScalar& c);
        ComplexScalar& operator=(ComplexScalar&& c) noexcept;
        ComplexScalar& operator=(const ScalarBase<ScalarType>& s);
        void operator<<=(int i) { real <<= i; imag<<= i; }
        void operator>>=(int i) { real >>= i; imag>>= i; }
        bool operator==(const ComplexScalar<ScalarType>& c) const;
        bool operator!=(const ComplexScalar<ScalarType>& c) const { return !(operator==(c)); }
        /* Helpers */
        void swap(ComplexScalar& c) noexcept;
        static inline ComplexScalar Zero();
        static inline ComplexScalar One();
        static inline ComplexScalar Two();
        static inline ComplexScalar Random();
        /* Getters */
        [[nodiscard]] ScalarType squaredNorm();
        [[nodiscard]] inline ScalarType norm();
        [[nodiscard]] ScalarType arg();
        [[nodiscard]] ComplexScalar conjugate() const noexcept { return ComplexScalar(real, -imag); }
        [[nodiscard]] const ScalarType& getReal() const { return real; }
        [[nodiscard]] const ScalarType& getImag() const { return imag; }
        [[nodiscard]] bool isZero() { return real.isZero() && imag.isZero(); }
        /* Setters */
        void setReal(const ScalarType& s) { real = s; }
        void setImag(const ScalarType& s) { imag = s; }
    };
    template<class ScalarType>
    std::ostream& operator<<(std::ostream& os, const ComplexScalar<ScalarType>& c);

    template<class ScalarType>
    ComplexScalar<ScalarType> operator+(
            const ComplexScalar<ScalarType>& c1, const ComplexScalar<ScalarType>& c2);

    template<class ScalarType>
    ComplexScalar<ScalarType> operator-(
            const ComplexScalar<ScalarType>& c1, const ComplexScalar<ScalarType>& c2);

    template<class ScalarType>
    ComplexScalar<ScalarType> operator*(
            const ComplexScalar<ScalarType>& c1, const ComplexScalar<ScalarType>& c2);

    template<class ScalarType>
    ComplexScalar<ScalarType> operator/(
            const ComplexScalar<ScalarType>& c1, const ComplexScalar<ScalarType>& c2);

    template<class ScalarType1, class ScalarType2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator+(
            const ComplexScalar<ScalarType1>& c, const ScalarBase<ScalarType2>& s);

    template<class ScalarType1, class ScalarType2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator-(
            const ComplexScalar<ScalarType1>& c, const ScalarBase<ScalarType2>& s);

    template<class ScalarType1, class ScalarType2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator*(
            const ComplexScalar<ScalarType1>& c, const ScalarBase<ScalarType2>& s);

    template<class ScalarType1, class ScalarType2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator/(
            const ComplexScalar<ScalarType1>& c, const ScalarBase<ScalarType2>& s);

    template<class ScalarType1, class ScalarType2>
    inline ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator+(
            const ScalarBase<ScalarType1>& s, const ComplexScalar<ScalarType2>& c);

    template<class ScalarType1, class ScalarType2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator-(
            const ScalarBase<ScalarType1>& s, const ComplexScalar<ScalarType2>& c);

    template<class ScalarType1, class ScalarType2>
    inline ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator*(
            const ScalarBase<ScalarType1>& s, const ComplexScalar<ScalarType2>& c);

    template<class ScalarType1, class ScalarType2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator/(
            const ScalarBase<ScalarType1>& s, const ComplexScalar<ScalarType2>& c);

    template<class ScalarType>
    ComplexScalar<ScalarType> operator<<(const ComplexScalar<ScalarType>& c, int i) {
        return ComplexScalar<ScalarType>(c.getReal() << i, c.getImag() << i);
    }

    template<class ScalarType>
    ComplexScalar<ScalarType> operator>>(const ComplexScalar<ScalarType>& c, int i) {
        return ComplexScalar<ScalarType>(c.getReal() >> i, c.getImag() >> i);
    }

    template<class ScalarType>
    void swap(ComplexScalar<ScalarType>& c1, ComplexScalar<ScalarType>& c2) noexcept { c1.swap(c2); }

    template<class ScalarType, class T>
    void operator+=(ComplexScalar<ScalarType>& c, const T& t) { c = c + t; }

    template<class ScalarType, class T>
    void operator-=(ComplexScalar<ScalarType>& c, const T& t) { c = c - t; }

    template<class ScalarType, class T>
    void operator*=(ComplexScalar<ScalarType>& c, const T& t) { c = c * t; }

    template<class ScalarType, class T>
    void operator/=(ComplexScalar<ScalarType>& c, const T& t) { c = c / t; }
}

#include "ComplexScalarImpl/ComplexScalarImpl.h"
