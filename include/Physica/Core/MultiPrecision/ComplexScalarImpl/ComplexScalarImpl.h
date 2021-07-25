/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.
 *
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
    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack>::ComplexScalar(Scalar<option, errorTrack> s1, Scalar<option, errorTrack> s2)
            : real(std::move(s1)), imag(std::move(s2)) {}

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack>::ComplexScalar(ComplexScalar&& c) noexcept
            : real(std::move(c.real)), imag(std::move(c.imag)) { Q_UNUSED(option) Q_UNUSED(errorTrack) }

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack>& ComplexScalar<option, errorTrack>::operator=(const ComplexScalar& c) {
        Q_UNUSED(option)
        Q_UNUSED(errorTrack)
        if(this == &c)
            return *this;
        real = c.real;
        imag = c.imag;
        return *this;
    }

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack>& ComplexScalar<option, errorTrack>::operator=(ComplexScalar&& c) noexcept {
        Q_UNUSED(option)
        Q_UNUSED(errorTrack)
        real = std::move(c.real);
        imag = std::move(c.imag);
        return *this;
    }

    template<ScalarOption option, bool errorTrack>
    bool ComplexScalar<option, errorTrack>::operator==(const ComplexScalar<option, errorTrack>& c) {
        return real == c.real && imag == c.imag;
    }

    template<ScalarOption option, bool errorTrack>
    bool ComplexScalar<option, errorTrack>::operator>(const ComplexScalar<option, errorTrack>& c) {
        return (square(real) + square(imag)) > (square(c.real) + square(c.imag));
    }

    template<ScalarOption option, bool errorTrack>
    bool ComplexScalar<option, errorTrack>::operator<(const ComplexScalar<option, errorTrack>& c) {
        return (square(real) + square(imag)) < (square(c.real) + square(c.imag));
    }

    template<ScalarOption option, bool errorTrack>
    void ComplexScalar<option, errorTrack>::swap(ComplexScalar& c) noexcept {
        Q_UNUSED(option)
        Q_UNUSED(errorTrack)
        Physica::Core::swap(real, c.real);
        Physica::Core::swap(imag, c.imag);
    }

    template<ScalarOption option, bool errorTrack>
    inline ComplexScalar<option, errorTrack> ComplexScalar<option, errorTrack>::getZero() {
        return ComplexScalar(Scalar<option, errorTrack>::getZero(), Scalar<option, errorTrack>::getZero());
    }

    template<ScalarOption option, bool errorTrack>
    inline ComplexScalar<option, errorTrack> ComplexScalar<option, errorTrack>::getOne() {
        return ComplexScalar(Scalar<option, errorTrack>::getOne(), Scalar<option, errorTrack>::getZero());
    }

    template<ScalarOption option, bool errorTrack>
    inline ComplexScalar<option, errorTrack> ComplexScalar<option, errorTrack>::getTwo() {
        return ComplexScalar(Scalar<option, errorTrack>::getTwo(), Scalar<option, errorTrack>::getZero());
    }

    template<ScalarOption option, bool errorTrack>
    inline ComplexScalar<option, errorTrack> ComplexScalar<option, errorTrack>::getRandom() {
        return ComplexScalar(randomScalar<option, errorTrack>(), randomScalar<option, errorTrack>());
    }

    template<ScalarOption option, bool errorTrack>
    inline Scalar<option, errorTrack> norm(const ComplexScalar<option, errorTrack>& c) {
        return sqrt(square(c.getReal()) + square(c.getImag()));
    }

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> arg(const ComplexScalar<option, errorTrack>& c) {
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

    template<ScalarOption option, bool errorTrack>
    std::ostream& operator<<(std::ostream& os, const ComplexScalar<option, errorTrack>& c) {
        const auto& imagine = c.getImag();
        return os << std::setprecision(10) << double(c.getReal())
                  << (imagine.isNegative() ? " - " : " + " ) << fabs(double(imagine)) << 'i' << std::setprecision(6);
    }

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> operator+(
            const ComplexScalar<option, errorTrack>& c1, const ComplexScalar<option, errorTrack>& c2) {
        return ComplexScalar<option, errorTrack>(c1.getReal() + c2.getReal(), c1.getImag() + c2.getImag());
    }

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> operator-(
            const ComplexScalar<option, errorTrack>& c1, const ComplexScalar<option, errorTrack>& c2) {
        return ComplexScalar<option, errorTrack>(c1.getReal() - c2.getReal(), c1.getImag() - c2.getImag());
    }

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> operator*(
            const ComplexScalar<option, errorTrack>& c1, const ComplexScalar<option, errorTrack>& c2) {
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
        return ComplexScalar<option, errorTrack>(ac - bd
                , (real_1 + real_2) * (imagine_1 + imagine_2) - ac - bd);
    }

    template<ScalarOption option, bool errorTrack>
    ComplexScalar<option, errorTrack> operator/(
            const ComplexScalar<option, errorTrack>& c1, const ComplexScalar<option, errorTrack>& c2) {
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
        return ComplexScalar<option, errorTrack>((ac + bd) / divisor
                , ((real_1 + imagine_1) * (real_2 - imagine_2) - ac + bd) / divisor);
    }

    template<ScalarOption option, bool errorTrack1, bool errorTrack2>
    ComplexScalar<option, errorTrack1 || errorTrack2> operator+(
            const ComplexScalar<option, errorTrack1>& c, const Scalar<option, errorTrack2>& s) {
        return ComplexScalar<option, errorTrack1 || errorTrack2>(c.getReal() + s, c.getImag());
    }

    template<ScalarOption option, bool errorTrack1, bool errorTrack2>
    ComplexScalar<option, errorTrack1 || errorTrack2> operator-(
            const ComplexScalar<option, errorTrack1>& c, const Scalar<option, errorTrack2>& s) {
        return ComplexScalar<option, errorTrack1 || errorTrack2>(c.getReal() - s, c.getImag());
    }

    template<ScalarOption option, bool errorTrack1, bool errorTrack2>
    ComplexScalar<option, errorTrack1 || errorTrack2> operator*(
            const ComplexScalar<option, errorTrack1>& c, const Scalar<option, errorTrack2>& s) {
        return ComplexScalar<option, errorTrack1 || errorTrack2>(c.getReal() * s, c.getImag() * s);
    }

    template<ScalarOption option, bool errorTrack1, bool errorTrack2>
    ComplexScalar<option, errorTrack1 || errorTrack2> operator/(
            const ComplexScalar<option, errorTrack1>& c, const Scalar<option, errorTrack2>& s) {
        return ComplexScalar<option, errorTrack1 || errorTrack2>(c.getReal() / s, c.getImag() / s);
    }

    template<ScalarOption option, bool errorTrack1, bool errorTrack2>
    ComplexScalar<option, errorTrack1 || errorTrack2> operator+(
            const Scalar<option, errorTrack1>& c, const ComplexScalar<option, errorTrack2>& s) {
        return ComplexScalar<option, errorTrack1 || errorTrack2>(c.getReal() + s, c.getImagine());
    }

    template<ScalarOption option, bool errorTrack1, bool errorTrack2>
    ComplexScalar<option, errorTrack1 || errorTrack2> operator-(
            const Scalar<option, errorTrack1>& c, const ComplexScalar<option, errorTrack2>& s) {
        return ComplexScalar<option, errorTrack1 || errorTrack2>(s - c.getReal(), c.getImagine());
    }

    template<ScalarOption option, bool errorTrack1, bool errorTrack2>
    ComplexScalar<option, errorTrack1 || errorTrack2> operator*(
            const Scalar<option, errorTrack1>& c, const ComplexScalar<option, errorTrack2>& s) {
        return ComplexScalar<option, errorTrack1 || errorTrack2>(c.getReal() * s, c.getImagine() * s);
    }

    template<ScalarOption option, bool errorTrack1, bool errorTrack2>
    ComplexScalar<option, errorTrack1 || errorTrack2> operator/(
            const Scalar<option, errorTrack1>& c, const ComplexScalar<option, errorTrack2>& s) {
        const auto& real = c.getReal();
        const auto& imagine = c.getImagine();
        const auto divisor = s * reciprocal(square(real) + square(imagine));
        return ComplexScalar<option, errorTrack1 || errorTrack2>(real * divisor, -imagine * divisor);
    }
}

#include "CElementaryFunction.h"

#endif