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
#ifndef PHYSICA_SCALARIMPL_H
#define PHYSICA_SCALARIMPL_H

#include <iomanip>
#include "ScalarArithmetic.h"
#include "ScalarAddSubExpression.h"
/**
 * This file is part of implementations of \Scalar.
 * Do not include this header file, include Scalar.h instead.
 */
namespace Physica::Core {
    //////////////////////////////////////////////Global//////////////////////////////////////////////
    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ScalarAddSubExpression<type, errorTrack1 || errorTrack2>
    operator+(const Scalar<type, errorTrack1>& s1, const Scalar<type, errorTrack2>& s2) {
        return ScalarAddSubExpression<type, errorTrack1 || errorTrack2>(s1, s2, false);
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ScalarAddSubExpression<type, errorTrack1 || errorTrack2>
    operator-(const Scalar<type, errorTrack1>& s1, const Scalar<type, errorTrack2>& s2) {
        return ScalarAddSubExpression<type, errorTrack1 || errorTrack2>(s1, s2, true);
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ScalarAddSubExpression<type, errorTrack1 || errorTrack2>
    operator+(const Scalar<type, errorTrack1>& s, ScalarAddSubExpression<type, errorTrack2>&& exp) {
        return ScalarAddSubExpression<type, errorTrack1 || errorTrack2>(s, std::move(exp), false);
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    ScalarAddSubExpression<type, errorTrack1 || errorTrack2>
    operator-(const Scalar<type, errorTrack1>& s, ScalarAddSubExpression<type, errorTrack2>&& exp) {
        return ScalarAddSubExpression<type, errorTrack1 || errorTrack2>(s, std::move(exp), true);
    }

    template<ScalarType type, bool errorTrack>
    std::ostream& operator<<(std::ostream& os, const Scalar<type, errorTrack>& s) {
        return os << std::setprecision(10) //10 is the max precision of double.
                  << double(s)
                  << std::setprecision(6); //6 is the default precision.
    }

    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> operator+(const Scalar<type, errorTrack>& s) {
        return Scalar<type, errorTrack>(s);
    }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    inline void operator+=(Scalar<type, errorTrack1>& s1
            , const Scalar<type, errorTrack2>& s2) { s1 = s1 + s2; }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    inline void operator-=(Scalar<type, errorTrack1>& s1
            , const Scalar<type, errorTrack2>& s2) { s1 = s1 - s2; }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    inline void operator*=(Scalar<type, errorTrack1>& s1
            , const Scalar<type, errorTrack2>& s2) { s1 = s1 * s2; }

    template<ScalarType type, bool errorTrack1, bool errorTrack2>
    inline void operator/=(Scalar<type, errorTrack1>& s1
            , const Scalar<type, errorTrack2>& s2) { s1 = s1 / s2; }

    template<ScalarType type, bool errorTrack>
    inline void operator^=(Scalar<type, errorTrack>& s1
            , const Scalar<type, errorTrack>& s2) { s1 = s1 ^ s2; }

    template<ScalarType type, bool errorTrack>
    inline void operator<<=(Scalar<type, errorTrack>& s
            , int bits) { s = s << bits; }

    template<ScalarType type, bool errorTrack>
    inline void operator>>=(Scalar<type, errorTrack>& s
            , int bits) { s = s >> bits; }
    /*!
     * The following two functions handle swap(s1, s2). Use swap of Scalars whose errorTrack is false by default.
     */
    template<ScalarType type, bool errorTrack>
    inline void swap(Scalar<type, false>& s1, Scalar<type, errorTrack>& s2) noexcept {
        s1.swap(s2);
    }

    template<ScalarType type, bool errorTrack>
    inline void swap(Scalar<type, true>& s1, Scalar<type, errorTrack>& s2) noexcept {
        s2.swap(s1);
    }
    ///////////////////////////////////////////MultiPrecision/////////////////////////////////////////
    /**
     * Returns true if s1 and s2 has the same sign. Both s1 and s2 do not equal to zero.
     * This function provide a quick sign check compare to using isPositive() and isNegative().
     */
    inline bool Scalar<MultiPrecision, false>::matchSign(
            const Scalar<MultiPrecision, false>& s1, const Scalar<MultiPrecision, false>& s2) {
        Q_ASSERT(!s1.isZero() && !s2.isZero());
        return (s1.length ^ s2.length) >= 0; //NOLINT Bitwise operator between two signed integer is intended.
    }

    inline Scalar<MultiPrecision, false>& operator++(Scalar<MultiPrecision, false>& s) {
        s += BasicConst::getInstance()._1;
        return s;
    }
    
    inline Scalar<MultiPrecision, false>& operator--(Scalar<MultiPrecision, false>& s) {
        s -= BasicConst::getInstance()._1;
        return s;
    }

    template<bool errorTrack>
    inline Scalar<MultiPrecision, errorTrack> operator++(Scalar<MultiPrecision, errorTrack>& s, int) { //NOLINT confusing-warning
        Scalar<MultiPrecision, errorTrack> temp(s);
        s += BasicConst::getInstance()._1;
        return temp;
    }

    template<bool errorTrack>
    inline Scalar<MultiPrecision, errorTrack> operator--(Scalar<MultiPrecision, errorTrack>& s, int) { //NOLINT confusing-warning
        Scalar<MultiPrecision, errorTrack> temp(s);
        s -= BasicConst::getInstance()._1;
        return temp;
    }
    //////////////////////////////////MultiPrecision-WithoutError///////////////////////////////////
    inline bool operator>=(const Scalar<MultiPrecision, false>& s1, const Scalar<MultiPrecision, false>& s2) {
        return !(s1 < s2);
    }

    inline bool operator<=(const Scalar<MultiPrecision, false>& s1, const Scalar<MultiPrecision, false>& s2) {
        return !(s1 > s2);
    }

    inline bool operator!= (const Scalar<MultiPrecision, false>& s1, const Scalar<MultiPrecision, false>& s2) {
        return !(s1 == s2);
    }
    ///////////////////////////////////MultiPrecision-WithError/////////////////////////////////////
    inline Scalar<MultiPrecision, false> Scalar<MultiPrecision, true>::getMaximum() const {
        return static_cast<const Scalar<MultiPrecision, false>&>(*this) + getAccuracy();
    }

    inline Scalar<MultiPrecision, false> Scalar<MultiPrecision, true>::getMinimum() const {
        return static_cast<const Scalar<MultiPrecision, false>&>(*this) - getAccuracy();
    }
    ///////////////////////////////////////////Float-Double////////////////////////////////////////////////
    /////////////////////////////////////////////Float////////////////////////////////////////////////
    template<bool errorTrack>
    inline Scalar<Float, errorTrack>& operator++(
            Scalar<Float, errorTrack>& s) {
        s += 1.0F;
        return s;
    }

    template<bool errorTrack>
    inline Scalar<Float, errorTrack>& operator--(
            Scalar<Float, errorTrack>& s) {
        s -= 1.0F;
        return s;
    }

    template<bool errorTrack>
    inline Scalar<Float, errorTrack> operator++( //NOLINT confusing-warning
            Scalar<Float, errorTrack>& s, int) {
        Scalar<Float, errorTrack> temp(s);
        s += 1.0F;
        return temp;
    }

    template<bool errorTrack>
    inline Scalar<Float, errorTrack> operator--( //NOLINT confusing-warning
            Scalar<Float, errorTrack>& s, int) {
        Scalar<Float, errorTrack> temp(s);
        s -= 1.0F;
        return temp;
    }
    ////////////////////////////////////////Float-WithoutError///////////////////////////////////////////
    inline Scalar<Float, false>::Scalar() : f(0) {}

    inline Scalar<Float, false>::Scalar(float f) : f(f) {}

    inline Scalar<Float, true> Scalar<Float, false>::operator*(const Scalar<Float, true>& s) const {
        return Scalar<Float, true>(f * s.f, f * s.getA());
    }

    inline Scalar<Float, true> Scalar<Float, false>::operator/(const Scalar<Float, true>& s) const {
        return Scalar<Float, true>(f / s.f, fabsf((s.f * s.getA()) / (s.f * (s.f - s.getA()))));
    }
    /*!
     * Set this scalar to its integer approximation which is closer to 0.
     * e.g. 5.6 -> 5, -4.8 -> -4, 0 -> 0.
     */
    inline void Scalar<Float, false>::toInteger() {
        modff(f, &f);
    }
    /////////////////////////////////////////Float-WithError////////////////////////////////////////////////
    inline Scalar<Float, true>::Scalar() : Scalar<Float, false>(), a(0) {}
    /*!
     * The abstract value of a equals to the accuracy.
     */
    inline Scalar<Float, true>::Scalar(float f, float a) : Scalar<Float, false>(f), a(fabsf(a)) {}

    inline Scalar<Float, true>::Scalar(const Scalar<Float, true>& s) : Scalar<Float, false>(s)  {
        a = s.a;
    }

    inline void Scalar<Float, true>::toInteger() {
        Scalar<Float, false>::toInteger();
        a = 0;
    }
    /////////////////////////////////////////////Double////////////////////////////////////////////////
    template<bool errorTrack>
    inline Scalar<Double, errorTrack>& operator++(Scalar<Double, errorTrack>& s) {
        s += 1.0;
        return s;
    }

    template<bool errorTrack>
    inline Scalar<Double, errorTrack>& operator--(Scalar<Double, errorTrack>& s) {
        s -= 1.0;
        return s;
    }

    template<bool errorTrack>
    inline Scalar<Double, errorTrack> operator++(Scalar<Double, errorTrack>& s, int) {  //NOLINT confusing-warning
        Scalar<Double, errorTrack> temp(s);
        s += 1.0;
        return temp;
    }

    template<bool errorTrack>
    inline Scalar<Double, errorTrack> operator--(Scalar<Double, errorTrack>& s, int) { //NOLINT confusing-warning
        Scalar<Double, errorTrack> temp(s);
        s -= 1.0;
        return temp;
    }
    ////////////////////////////////////////Double-WithoutError///////////////////////////////////////////
    inline Scalar<Double, false>::Scalar() : d(0) {}

    inline Scalar<Double, false>::Scalar(double d) : d(d) {}

    inline Scalar<Double, true> Scalar<Double, false>::operator*(const Scalar<Double, true>& s) const {
        return Scalar<Double, true>(d * s.d, d * s.getA());
    }

    inline Scalar<Double, true> Scalar<Double, false>::operator/(const Scalar<Double, true>& s) const {
        return Scalar<Double, true>(d / s.d, fabs((s.d * s.getA()) / (s.d * (s.d - s.getA()))));
    }
    /*!
     * Set this scalar to its integer approximation which is closer to 0.
     * e.g. 5.6 -> 5, -4.8 -> -4, 0 -> 0.
     */
    inline void Scalar<Double, false>::toInteger() {
        modf(d, &d);
    }
    /////////////////////////////////////////Double-WithError///////////////////////////////////////////
    inline Scalar<Double, true>::Scalar() : Scalar<Double, false>(), a(0) {}
    /*!
     * The abstract value of a equals to the accuracy.
     */
    inline Scalar<Double, true>::Scalar(double d, double a) : Scalar<Double, false>(d), a(fabs(a)) {}

    inline Scalar<Double, true>::Scalar(const Scalar<Double, true>& s) : Scalar<Double, false>(s) {
        a = s.a;
    }

    inline void Scalar<Double, true>::toInteger() {
        Scalar<Double, false>::toInteger();
        a = 0;
    }
}

#include "ScalarArithmetic.h"
#include "ElementaryFunction.h"
#include "Operation/Pow.h"

namespace Physica::Core {
    template<ScalarType type>
    inline Scalar<type, false> operator^(
            const Scalar<type, false>& s1, const Scalar<type, false>& s2) {
        return Scalar<type, false>(std::pow(s1.getTrivial(), s2.getTrivial()));
    }

    template<ScalarType type>
    inline Scalar<type, true> operator^(
            const Scalar<type, true>& s1, const Scalar<type, false>& s2) {
        const auto result = std::pow(s1.getTrivial(), s2.getTrivial());
        return Scalar<type, true>(result, std::pow(s1.getTrivial(), s2.getTrivial() + s2.getA()) - result);
    }

    template<ScalarType type>
    inline Scalar<type, true> operator^(
            const Scalar<type, false>& s1, const Scalar<type, true>& s2) {
        const auto result = std::pow(s1.getTrivial(), s2.getTrivial());
        return Scalar<type, true>(result
                , std::pow(s1.getTrivial() + (s1.getTrivial() > 1 ? s1.getA() : -s1.getA()), s2.getTrivial()) - result);
    }

    template<ScalarType type>
    inline Scalar<type, true> operator^(
            const Scalar<type, true>& s1, const Scalar<type, true>& s2) {
        const auto result = std::pow(s1.getTrivial(), s2.getTrivial());
        return Scalar<type, true>(result
                , std::pow(s1.getTrivial() + (s1.getTrivial() > 1 ? s1.getA() : -s1.getA())
                        , s2.getTrivial() + s2.getA()) - result);
    }

    template<>
    inline Scalar<MultiPrecision, false> operator^(
            const Scalar<MultiPrecision, false>& s1, const Scalar<MultiPrecision, false>& s2) {
        return s1.isInteger() ? powScalar(s1, s2) : exp(ln(s1) * s2);
    }

    template<>
    inline Scalar<MultiPrecision, true> operator^(
            const Scalar<MultiPrecision, true>& s1, const Scalar<MultiPrecision, false>& s2) {
        return exp(ln(s1) * s2);
    }

    template<>
    inline Scalar<MultiPrecision, true> operator^(
            const Scalar<MultiPrecision, false>& s1, const Scalar<MultiPrecision, true>& s2) {
        return exp(ln(s1) * s2);
    }

    template<>
    inline Scalar<MultiPrecision, true> operator^(
            const Scalar<MultiPrecision, true>& s1, const Scalar<MultiPrecision, true>& s2) {
        return exp(ln(s1) * s2);
    }
}

#endif