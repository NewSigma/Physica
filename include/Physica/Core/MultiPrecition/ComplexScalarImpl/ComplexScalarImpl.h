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
#ifndef PHYSICA_COMPLEXSCALARIMPL_H
#define PHYSICA_COMPLEXSCALARIMPL_H
/*!
 * This file is part of implementations of \ComplexScalar.
 * Do not include this header file, include ComplexScalar.h instead.
 */
namespace Physica::Core {
    //!Optimize: maybe use && to avoid unnecessary move.
    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack>::ComplexScalar(Scalar<type, errorTrack> s1, Scalar<type, errorTrack> s2)
            : real(std::move(s1)), imag(std::move(s2)) {}

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack>::ComplexScalar(ComplexScalar&& c) noexcept
            : real(std::move(c.real)), imag(std::move(c.imag)) { Q_UNUSED(type) Q_UNUSED(errorTrack) }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack>& ComplexScalar<type, errorTrack>::operator=(const ComplexScalar& c) {
        Q_UNUSED(type)
        Q_UNUSED(errorTrack)
        if(this == &c)
            return *this;
        real = c.real;
        imag = c.imag;
        return *this;
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack>& ComplexScalar<type, errorTrack>::operator=(ComplexScalar&& c) noexcept {
        Q_UNUSED(type)
        Q_UNUSED(errorTrack)
        real = std::move(c.real);
        imag = std::move(c.imag);
        return *this;
    }

    template<ScalarType type, bool errorTrack>
    bool ComplexScalar<type, errorTrack>::operator==(const ComplexScalar<type, errorTrack>& c) {
        return real == c.real && imag == c.imag;
    }

    template<ScalarType type, bool errorTrack>
    bool ComplexScalar<type, errorTrack>::operator>(const ComplexScalar<type, errorTrack>& c) {
        return (square(real) + square(imag)) > (square(c.real) + square(c.imag));
    }

    template<ScalarType type, bool errorTrack>
    bool ComplexScalar<type, errorTrack>::operator<(const ComplexScalar<type, errorTrack>& c) {
        return (square(real) + square(imag)) < (square(c.real) + square(c.imag));
    }

    template<ScalarType type, bool errorTrack>
    void ComplexScalar<type, errorTrack>::swap(ComplexScalar& c) noexcept {
        Q_UNUSED(type)
        Q_UNUSED(errorTrack)
        Physica::Core::swap(real, c.real);
        Physica::Core::swap(imag, c.imag);
    }

    template<ScalarType type, bool errorTrack>
    inline ComplexScalar<type, errorTrack> ComplexScalar<type, errorTrack>::getZero() {
        return ComplexScalar(Scalar<type, errorTrack>::getZero(), Scalar<type, errorTrack>::getZero());
    }

    template<ScalarType type, bool errorTrack>
    inline ComplexScalar<type, errorTrack> ComplexScalar<type, errorTrack>::getOne() {
        return ComplexScalar(Scalar<type, errorTrack>::getOne(), Scalar<type, errorTrack>::getZero());
    }

    template<ScalarType type, bool errorTrack>
    inline ComplexScalar<type, errorTrack> ComplexScalar<type, errorTrack>::getTwo() {
        return ComplexScalar(Scalar<type, errorTrack>::getTwo(), Scalar<type, errorTrack>::getZero());
    }

    template<ScalarType type, bool errorTrack>
    inline ComplexScalar<type, errorTrack> ComplexScalar<type, errorTrack>::getRandom() {
        return ComplexScalar(randomScalar<type, errorTrack>(), randomScalar<type, errorTrack>());
    }

    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> norm(const ComplexScalar<type, errorTrack>& c) {
        return sqrt(square(c.getReal()) + square(c.getImag()));
    }

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> arg(const ComplexScalar<type, errorTrack>& c) {
        const auto& real = c.getReal();
        const auto& imagine  = c.getImag();
        auto result = arctan(imagine / real);
        //arctan is defined on [-Pi / 2, Pi / 2], judging the quadrant is necessary.
        if(real.isPositive()) {
            if(imagine.isNegative())
                result += MathConst::getInstance().PI << 1;
        }
        else
            result += MathConst::getInstance().PI;
        return result;
    }

    template<ScalarType type, bool errorTrack>
    std::ostream& operator<<(std::ostream& os, const ComplexScalar<type, errorTrack>& c) {
        const auto& imagine = c.getImag();
        return os << std::setprecision(10) << double(c.getReal())
                  << (imagine.isNegative() ? " - " : "+" )<< 'i' << fabs(double(imagine)) << std::setprecision(6);
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> operator+(
            const ComplexScalar<type, errorTrack>& c1, const ComplexScalar<type, errorTrack>& c2) {
        return ComplexScalar<type, errorTrack>(c1.getReal() + c2.getReal(), c1.getImag() + c2.getImag());
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> operator-(
            const ComplexScalar<type, errorTrack>& c1, const ComplexScalar<type, errorTrack>& c2) {
        return ComplexScalar<type, errorTrack>(c1.getReal() - c2.getReal(), c1.getImag() - c2.getImag());
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> operator*(
            const ComplexScalar<type, errorTrack>& c1, const ComplexScalar<type, errorTrack>& c2) {
        const auto& real_1 = c1.getReal();
        const auto& imagine_1 = c1.getImag();
        const auto& real_2 = c2.getReal();
        const auto& imagine_2 = c2.getImag();
        /*
         * Optimize:
         * Use (a + ib)(c + id) = (ac - bd) + i((a + b)(c + d) - ac - bd)
         * instead of (a + ib)(c + id) = (ac - bd) + i(ad + bc) to avoid multiply.
         * But it is unclear if this method is useful to every machine.
         * May be add checks and use Config.h to determine which method to use.
         */
        const auto ac = real_1 * real_2;
        const auto bd = imagine_1 * imagine_2;
        return ComplexScalar<type, errorTrack>(ac - bd
                , (real_1 + real_2) * (imagine_1 + imagine_2) - ac - bd);
    }

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> operator/(
            const ComplexScalar<type, errorTrack>& c1, const ComplexScalar<type, errorTrack>& c2) {
        const auto& real_1 = c1.getReal();
        const auto& imagine_1 = c1.getImag();
        const auto& real_2 = c2.getReal();
        const auto& imagine_2 = c2.getImag();
        /*
         * Optimize: Using the same method with operator*().
         */
        const auto ac = real_1 * real_2;
        const auto bd = imagine_1 * imagine_2;
        const auto divisor = square(real_2) + square(imagine_2);
        return ComplexScalar<type, errorTrack>((ac + bd) / divisor
                , ((real_1 + imagine_1) * (real_2 - imagine_2) - ac + bd) / divisor);
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator+(
            const ComplexScalar<type, errorTrack1>& c, const Scalar<type, errorTrack2>& s) {
        return ComplexScalar<type, errorTrack1 | errorTrack2>(c.getReal() + s, c.getImag());
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator-(
            const ComplexScalar<type, errorTrack1>& c, const Scalar<type, errorTrack2>& s) {
        return ComplexScalar<type, errorTrack1 | errorTrack2>(c.getReal() - s, c.getImag());
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator*(
            const ComplexScalar<type, errorTrack1>& c, const Scalar<type, errorTrack2>& s) {
        return ComplexScalar<type, errorTrack1 | errorTrack2>(c.getReal() * s, c.getImag() * s);
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator/(
            const ComplexScalar<type, errorTrack1>& c, const Scalar<type, errorTrack2>& s) {
        return ComplexScalar<type, errorTrack1 | errorTrack2>(c.getReal() / s, c.getImag() / s);
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator+(
            const Scalar<type, errorTrack1>& c, const ComplexScalar<type, errorTrack2>& s) {
        return ComplexScalar<type, errorTrack1 | errorTrack2>(c.getReal() + s, c.getImagine());
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator-(
            const Scalar<type, errorTrack1>& c, const ComplexScalar<type, errorTrack2>& s) {
        return ComplexScalar<type, errorTrack1 | errorTrack2>(s - c.getReal(), c.getImagine());
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator*(
            const Scalar<type, errorTrack1>& c, const ComplexScalar<type, errorTrack2>& s) {
        return ComplexScalar<type, errorTrack1 | errorTrack2>(c.getReal() * s, c.getImagine() * s);
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ComplexScalar<type, errorTrack1 | errorTrack2> operator/(
            const Scalar<type, errorTrack1>& c, const ComplexScalar<type, errorTrack2>& s) {
        const auto& real = c.getReal();
        const auto& imagine = c.getImagine();
        const auto divisor = s * reciprocal(square(real) + square(imagine));
        return ComplexScalar<type, errorTrack1 | errorTrack2>(real * divisor, -imagine * divisor);
    }
}

#include "CElementaryFunction.h"

#endif