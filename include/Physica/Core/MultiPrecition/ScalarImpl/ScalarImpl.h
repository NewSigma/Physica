/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_SCALARIMPL_H
#define PHYSICA_SCALARIMPL_H

#include <iomanip>
#include "Util/Bitwise.h"
/*!
 * This file is part of implementations of \Scalar.
 * Do not include this header file, include Scalar.h instead.
 */
namespace Physica::Core {
    //////////////////////////////////////////////Global//////////////////////////////////////////////
    template<ScalarType type, bool errorTrack>
    std::ostream& operator<<(std::ostream& os, const Scalar<type, errorTrack>& n) {
        return os << std::setprecision(10) //10 is the max precision of double.
                  << double(n)
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
    //!Reference: MaTHmu Project Group.计算机代数系统的数学原理[M].Beijing: TsingHua University Press, 2009.45
    template<bool errorTrack1, bool errorTrack2>
    inline Scalar<MultiPrecision, errorTrack1 | errorTrack2> operator^(
            const Scalar<MultiPrecision, errorTrack1>& s1,
            const Scalar<MultiPrecision, errorTrack2>& s2) {
        return s1.isInteger() ? powScalar(s1, s2) : exp(ln(s1) * s2);
    }

    inline Scalar<MultiPrecision, false>& operator++(Scalar<MultiPrecision, false>& s) {
        s += BasicConst::getInstance().get_1();
        return s;
    }
    
    inline Scalar<MultiPrecision, false>& operator--(Scalar<MultiPrecision, false>& s) {
        s -= BasicConst::getInstance().get_1();
        return s;
    }

    template<bool errorTrack>
    inline Scalar<MultiPrecision, errorTrack> operator++(Scalar<MultiPrecision, errorTrack>& s, int) { //NOLINT confusing-warning
        Scalar<MultiPrecision, errorTrack> temp(s);
        s += BasicConst::getInstance().get_1();
        return temp;
    }

    template<bool errorTrack>
    inline Scalar<MultiPrecision, errorTrack> operator--(Scalar<MultiPrecision, errorTrack>& s, int) { //NOLINT confusing-warning
        Scalar<MultiPrecision, errorTrack> temp(s);
        s -= BasicConst::getInstance().get_1();
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

    inline Scalar<Float, true> Scalar<Float, false>::operator+(const Scalar<Float, true>& s) const {
        return Scalar<Float, true>(f + s.f, s.getA());
    }

    inline Scalar<Float, true> Scalar<Float, false>::operator-(const Scalar<Float, true>& s) const {
        return Scalar<Float, true>(f - s.f, s.getA());
    }

    inline Scalar<Float, true> Scalar<Float, false>::operator*(const Scalar<Float, true>& s) const {
        return Scalar<Float, true>(f * s.f, f * s.getA());
    }

    inline Scalar<Float, true> Scalar<Float, false>::operator/(const Scalar<Float, true>& s) const {
        return Scalar<Float, true>(f / s.f, fabsf((s.f * s.getA()) / (s.f * (s.f - s.getA()))));
    }
    /////////////////////////////////////////Float-WithError////////////////////////////////////////////////
    inline Scalar<Float, true>::Scalar() : Scalar<Float, false>(), a(0) {}
    /*!
     * The abstract value of a equals to the accuracy.
     */
    inline Scalar<Float, true>::Scalar(float f, float a) : Scalar<Float, false>(f), a(fabsf(a)) {}
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

    inline Scalar<Double, true> Scalar<Double, false>::operator+(const Scalar<Double, true>& s) const {
        return Scalar<Double, true>(d + s.d, s.getA());
    }

    inline Scalar<Double, true> Scalar<Double, false>::operator-(const Scalar<Double, true>& s) const {
        return Scalar<Double, true>(d - s.d, s.getA());
    }

    inline Scalar<Double, true> Scalar<Double, false>::operator*(const Scalar<Double, true>& s) const {
        return Scalar<Double, true>(d * s.d, d * s.getA());
    }

    inline Scalar<Double, true> Scalar<Double, false>::operator/(const Scalar<Double, true>& s) const {
        return Scalar<Double, true>(d / s.d, fabs((s.d * s.getA()) / (s.d * (s.d - s.getA()))));
    }
    /////////////////////////////////////////Double-WithError///////////////////////////////////////////
    inline Scalar<Double, true>::Scalar() : Scalar<Double, false>(), a(0) {}
    /*!
     * The abstract value of a equals to the accuracy.
     */
    inline Scalar<Double, true>::Scalar(double d, double a) : Scalar<Double, false>(d), a(fabs(a)) {}
}

#endif