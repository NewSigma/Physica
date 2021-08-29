/*
 * Copyright 2020-2021 WeiBo He.
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
#pragma once

namespace Physica::Core {
    template<class AnyScalar>
    ComplexScalar<AnyScalar>::ComplexScalar(const ScalarBase<AnyScalar>& real_)
            : real(real_.getDerived()), imag(AnyScalar::Zero()) {}

    template<class AnyScalar>
    ComplexScalar<AnyScalar>::ComplexScalar(const ScalarBase<AnyScalar>& real_, const ScalarBase<AnyScalar>& imag_)
            : real(real_.getDerived()), imag(imag_.getDerived()) {}

    template<class AnyScalar>
    ComplexScalar<AnyScalar>::ComplexScalar(std::initializer_list<AnyScalar> list) {
        assert(list.size() == 2);
        auto ite = list.begin();
        real = *ite;
        ++ite;
        imag = *ite;
    }

    template<class AnyScalar>
    ComplexScalar<AnyScalar>::ComplexScalar(ComplexScalar&& c) noexcept
            : real(std::move(c.real)), imag(std::move(c.imag)) {}

    template<class AnyScalar>
    ComplexScalar<AnyScalar>& ComplexScalar<AnyScalar>::operator=(const ComplexScalar& c) {
        if(this == &c)
            return *this;
        real = c.real;
        imag = c.imag;
        return *this;
    }

    template<class AnyScalar>
    ComplexScalar<AnyScalar>& ComplexScalar<AnyScalar>::operator=(ComplexScalar&& c) noexcept {
        real = std::move(c.real);
        imag = std::move(c.imag);
        return *this;
    }

    template<class AnyScalar>
    ComplexScalar<AnyScalar>& ComplexScalar<AnyScalar>::operator=(const ScalarBase<AnyScalar>& s) {
        real = s.getDerived();
        imag = AnyScalar::Zero();
        return *this;
    }

    template<class AnyScalar>
    bool ComplexScalar<AnyScalar>::operator==(const ComplexScalar<AnyScalar>& c) {
        return real == c.real && imag == c.imag;
    }

    template<class AnyScalar>
    bool ComplexScalar<AnyScalar>::operator>(const ComplexScalar<AnyScalar>& c) {
        return (square(real) + square(imag)) > (square(c.real) + square(c.imag));
    }

    template<class AnyScalar>
    bool ComplexScalar<AnyScalar>::operator<(const ComplexScalar<AnyScalar>& c) {
        return (square(real) + square(imag)) < (square(c.real) + square(c.imag));
    }

    template<class AnyScalar>
    void ComplexScalar<AnyScalar>::swap(ComplexScalar& c) noexcept {
        Physica::Core::swap(real, c.real);
        Physica::Core::swap(imag, c.imag);
    }

    template<class AnyScalar>
    inline ComplexScalar<AnyScalar> ComplexScalar<AnyScalar>::Zero() {
        return ComplexScalar(AnyScalar::Zero(), AnyScalar::Zero());
    }

    template<class AnyScalar>
    inline ComplexScalar<AnyScalar> ComplexScalar<AnyScalar>::One() {
        return ComplexScalar(AnyScalar::One(), AnyScalar::Zero());
    }

    template<class AnyScalar>
    inline ComplexScalar<AnyScalar> ComplexScalar<AnyScalar>::Two() {
        return ComplexScalar(AnyScalar::Two(), AnyScalar::Zero());
    }

    template<class AnyScalar>
    inline ComplexScalar<AnyScalar> ComplexScalar<AnyScalar>::Random() {
        return ComplexScalar(randomScalar<AnyScalar>(), randomScalar<AnyScalar>());
    }

    template<class AnyScalar>
    inline AnyScalar ComplexScalar<AnyScalar>::norm() {
        return sqrt(square(real) + square(imag));
    }

    template<class AnyScalar>
    AnyScalar ComplexScalar<AnyScalar>::arg() {
        auto result = arctan(imag / real);
        //arctan is defined on [-Pi / 2, Pi / 2], judging the quadrant is necessary.
        if(real.isPositive()) {
            if(imag.isNegative())
                result += MathConst::getInstance().PI << 1;
        }
        else
            result += MathConst::getInstance().PI;
        return result;
    }

    template<class AnyScalar>
    std::ostream& operator<<(std::ostream& os, const ComplexScalar<AnyScalar>& c) {
        const auto& imagine = c.getImag();
        return os << std::setprecision(10) << double(c.getReal())
                  << (imagine.isNegative() ? " - " : " + " ) << fabs(double(imagine)) << 'i' << std::setprecision(6);
    }

    template<class AnyScalar>
    ComplexScalar<AnyScalar> operator+(
            const ComplexScalar<AnyScalar>& c1, const ComplexScalar<AnyScalar>& c2) {
        return ComplexScalar<AnyScalar>(c1.getReal() + c2.getReal(), c1.getImag() + c2.getImag());
    }

    template<class AnyScalar>
    ComplexScalar<AnyScalar> operator-(
            const ComplexScalar<AnyScalar>& c1, const ComplexScalar<AnyScalar>& c2) {
        return ComplexScalar<AnyScalar>(c1.getReal() - c2.getReal(), c1.getImag() - c2.getImag());
    }

    template<class AnyScalar>
    ComplexScalar<AnyScalar> operator*(
            const ComplexScalar<AnyScalar>& c1, const ComplexScalar<AnyScalar>& c2) {
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
        return ComplexScalar<AnyScalar>(ac - bd
                , (real_1 + real_2) * (imagine_1 + imagine_2) - ac - bd);
    }

    template<class AnyScalar>
    ComplexScalar<AnyScalar> operator/(
            const ComplexScalar<AnyScalar>& c1, const ComplexScalar<AnyScalar>& c2) {
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
        return ComplexScalar<AnyScalar>((ac + bd) / divisor
                , ((real_1 + imagine_1) * (real_2 - imagine_2) - ac + bd) / divisor);
    }

    template<class AnyScalar1, class AnyScalar2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type> operator+(
            const ComplexScalar<AnyScalar1>& c, const ScalarBase<AnyScalar2>& s) {
        return ComplexScalar<typename Internal::BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type>(c.getReal() + s.getDerived(), c.getImag());
    }

    template<class AnyScalar1, class AnyScalar2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type> operator-(
            const ComplexScalar<AnyScalar1>& c, const ScalarBase<AnyScalar2>& s) {
        return ComplexScalar<typename Internal::BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type>(c.getReal() - s.getDerived(), c.getImag());
    }

    template<class AnyScalar1, class AnyScalar2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type> operator*(
            const ComplexScalar<AnyScalar1>& c, const ScalarBase<AnyScalar2>& s) {
        return ComplexScalar<typename Internal::BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type>(c.getReal() * s.getDerived(), c.getImag() * s.getDerived());
    }

    template<class AnyScalar1, class AnyScalar2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type> operator/(
            const ComplexScalar<AnyScalar1>& c, const ScalarBase<AnyScalar2>& s) {
        return ComplexScalar<typename Internal::BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type>(c.getReal() / s.getDerived(), c.getImag() / s.getDerived());
    }

    template<class AnyScalar1, class AnyScalar2>
    inline ComplexScalar<typename Internal::BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type> operator+(
            const ScalarBase<AnyScalar1>& s, const ComplexScalar<AnyScalar2>& c) {
        return c + s;
    }

    template<class AnyScalar1, class AnyScalar2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type> operator-(
            const ScalarBase<AnyScalar1>& s, const ComplexScalar<AnyScalar2>& c) {
        return ComplexScalar<typename Internal::BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type>(s.getDerived() - c.getReal(), c.getImagine());
    }

    template<class AnyScalar1, class AnyScalar2>
    inline ComplexScalar<typename Internal::BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type> operator*(
            const ScalarBase<AnyScalar1>& s, const ComplexScalar<AnyScalar2>& c) {
        return c * s;
    }

    template<class AnyScalar1, class AnyScalar2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type> operator/(
            const ScalarBase<AnyScalar1>& s, const ComplexScalar<AnyScalar2>& c) {
        const auto& real = c.getReal();
        const auto& imagine = c.getImagine();
        const auto divisor = s * reciprocal(square(real) + square(imagine));
        return ComplexScalar<typename Internal::BinaryScalarOpReturnType<AnyScalar1, AnyScalar2>::Type>(real * divisor, -imagine * divisor);
    }
}

#include "CElementaryFunction.h"
