/*
 * Copyright 2020-2022 WeiBo He.
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
    template<class ScalarType>
    ComplexScalar<ScalarType>::ComplexScalar(const ScalarType& real_)
            : real(real_), imag(ScalarType::Zero()) {}

    template<class ScalarType>
    ComplexScalar<ScalarType>::ComplexScalar(const ScalarType& real_, const ScalarType& imag_)
            : real(real_), imag(imag_) {}

    template<class ScalarType>
    ComplexScalar<ScalarType>::ComplexScalar(std::initializer_list<ScalarType> list) {
        assert(list.size() == 2);
        auto ite = list.begin();
        real = *ite;
        ++ite;
        imag = *ite;
    }

    template<class ScalarType>
    ComplexScalar<ScalarType>::ComplexScalar(ComplexScalar&& c) noexcept
            : real(std::move(c.real)), imag(std::move(c.imag)) {}

    template<class ScalarType>
    ComplexScalar<ScalarType>& ComplexScalar<ScalarType>::operator=(const ComplexScalar& c) {
        if(this == &c)
            return *this;
        real = c.real;
        imag = c.imag;
        return *this;
    }

    template<class ScalarType>
    ComplexScalar<ScalarType>& ComplexScalar<ScalarType>::operator=(ComplexScalar&& c) noexcept {
        real = std::move(c.real);
        imag = std::move(c.imag);
        return *this;
    }

    template<class ScalarType>
    ComplexScalar<ScalarType>& ComplexScalar<ScalarType>::operator=(const ScalarBase<ScalarType>& s) {
        real = s.getDerived();
        imag = ScalarType::Zero();
        return *this;
    }

    template<class ScalarType>
    bool ComplexScalar<ScalarType>::operator==(const ComplexScalar<ScalarType>& c) const {
        return real == c.real && imag == c.imag;
    }

    template<class ScalarType>
    void ComplexScalar<ScalarType>::swap(ComplexScalar& c) noexcept {
        Physica::Core::swap(real, c.real);
        Physica::Core::swap(imag, c.imag);
    }

    template<class ScalarType>
    inline ComplexScalar<ScalarType> ComplexScalar<ScalarType>::Zero() {
        return ComplexScalar(ScalarType::Zero(), ScalarType::Zero());
    }

    template<class ScalarType>
    inline ComplexScalar<ScalarType> ComplexScalar<ScalarType>::One() {
        return ComplexScalar(ScalarType::One(), ScalarType::Zero());
    }

    template<class ScalarType>
    inline ComplexScalar<ScalarType> ComplexScalar<ScalarType>::Two() {
        return ComplexScalar(ScalarType::Two(), ScalarType::Zero());
    }

    template<class ScalarType>
    inline ComplexScalar<ScalarType> ComplexScalar<ScalarType>::Random() {
        return ComplexScalar(randomScalar<ScalarType>(), randomScalar<ScalarType>());
    }

    template<class ScalarType>
    ScalarType ComplexScalar<ScalarType>::squaredNorm() const {
        return square(real) + square(imag);
    }

    template<class ScalarType>
    inline ScalarType ComplexScalar<ScalarType>::norm() const {
        return sqrt(squaredNorm());
    }

    template<class ScalarType>
    ScalarType ComplexScalar<ScalarType>::arg() const {
        using RealType = typename ScalarType::RealType;
        auto result = arctan(imag / real);
        //arctan is defined on [-Pi / 2, Pi / 2], judging the quadrant is necessary.
        if(real.isPositive()) {
            if(imag.isNegative())
                result += RealType(M_PI_2);
        }
        else
            result += RealType(M_PI);
        return result;
    }

    template<class ScalarType>
    std::ostream& operator<<(std::ostream& os, const ComplexScalar<ScalarType>& c) {
        const auto& imagine = c.getImag();
        return os << std::setprecision(10) << double(c.getReal())
                  << (imagine.isNegative() ? " - " : " + " ) << fabs(double(imagine)) << 'i' << std::setprecision(6);
    }

    template<class ScalarType>
    ComplexScalar<ScalarType> operator+(
            const ComplexScalar<ScalarType>& c1, const ComplexScalar<ScalarType>& c2) {
        return ComplexScalar<ScalarType>(c1.getReal() + c2.getReal(), c1.getImag() + c2.getImag());
    }

    template<class ScalarType>
    ComplexScalar<ScalarType> operator-(
            const ComplexScalar<ScalarType>& c1, const ComplexScalar<ScalarType>& c2) {
        return ComplexScalar<ScalarType>(c1.getReal() - c2.getReal(), c1.getImag() - c2.getImag());
    }

    template<class ScalarType>
    ComplexScalar<ScalarType> operator*(
            const ComplexScalar<ScalarType>& c1, const ComplexScalar<ScalarType>& c2) {
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
        return ComplexScalar<ScalarType>(ac - bd
                , (real_1 + imagine_1) * (real_2 + imagine_2) - ac - bd);
    }

    template<class ScalarType>
    ComplexScalar<ScalarType> operator/(
            const ComplexScalar<ScalarType>& c1, const ComplexScalar<ScalarType>& c2) {
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
        return ComplexScalar<ScalarType>((ac + bd) / divisor
                , ((real_1 + imagine_1) * (real_2 - imagine_2) - ac + bd) / divisor);
    }

    template<class ScalarType1, class ScalarType2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator+(
            const ComplexScalar<ScalarType1>& c, const ScalarBase<ScalarType2>& s) {
        return ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type>(c.getReal() + s.getDerived(), c.getImag());
    }

    template<class ScalarType1, class ScalarType2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator-(
            const ComplexScalar<ScalarType1>& c, const ScalarBase<ScalarType2>& s) {
        return ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type>(c.getReal() - s.getDerived(), c.getImag());
    }

    template<class ScalarType1, class ScalarType2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator*(
            const ComplexScalar<ScalarType1>& c, const ScalarBase<ScalarType2>& s) {
        return ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type>(c.getReal() * s.getDerived(), c.getImag() * s.getDerived());
    }

    template<class ScalarType1, class ScalarType2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator/(
            const ComplexScalar<ScalarType1>& c, const ScalarBase<ScalarType2>& s) {
        return ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type>(c.getReal() / s.getDerived(), c.getImag() / s.getDerived());
    }

    template<class ScalarType1, class ScalarType2>
    inline ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator+(
            const ScalarBase<ScalarType1>& s, const ComplexScalar<ScalarType2>& c) {
        return c + s;
    }

    template<class ScalarType1, class ScalarType2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator-(
            const ScalarBase<ScalarType1>& s, const ComplexScalar<ScalarType2>& c) {
        return ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type>(s.getDerived() - c.getReal(), c.getImag());
    }

    template<class ScalarType1, class ScalarType2>
    inline ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator*(
            const ScalarBase<ScalarType1>& s, const ComplexScalar<ScalarType2>& c) {
        return c * s;
    }

    template<class ScalarType1, class ScalarType2>
    ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type> operator/(
            const ScalarBase<ScalarType1>& s, const ComplexScalar<ScalarType2>& c) {
        const auto& real = c.getReal();
        const auto& imagine = c.getImagine();
        const auto divisor = s * reciprocal(square(real) + square(imagine));
        return ComplexScalar<typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type>(real * divisor, -imagine * divisor);
    }
}

#include "CElementaryFunction.h"
