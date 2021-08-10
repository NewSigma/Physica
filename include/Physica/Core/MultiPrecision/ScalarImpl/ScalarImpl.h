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

#include <iomanip>
#include "ScalarArithmetic.h"
/**
 * This file is part of implementations of \Scalar.
 * Do not include this header file, include \file Scalar.h instead.
 */
namespace Physica::Core {
    namespace Internal {
        /**
         * Returns true if s1 and s2 has the same sign. Both s1 and s2 do not equal to zero.
         * This function provide a quick sign check compare to using isPositive() and isNegative().
         */
        inline bool AbstractScalar<MultiPrecision>::matchSign(
                const AbstractScalar<MultiPrecision>& s1, const AbstractScalar<MultiPrecision>& s2) {
            Q_ASSERT(!s1.isZero() && !s2.isZero());
            return (s1.length ^ s2.length) >= 0; //NOLINT Bitwise operator between two signed integer is intended.
        }
        /**
         * Cut zeros from the beginning.
         */
        inline void AbstractScalar<MultiPrecision>::cutZero() {
            const int size = getSize();
            int id = size - 1;
            while(byte[id] == 0 && id > 0)
                --id;
            ++id;

            if(id != size) {
                int shorten = size - id;
                byte = reinterpret_cast<MPUnit*>(realloc(byte, id * sizeof(MPUnit)));
                length = length > 0 ? id : -id;
                auto temp = power;
                power = byte[id - 1] != 0 ? (temp - shorten) : 0;
            }
        }
    }
    //////////////////////////////////////////////Global//////////////////////////////////////////////
    template<ScalarOption option, bool errorTrack>
    std::ostream& operator<<(std::ostream& os, const Scalar<option, errorTrack>& s) {
        return os << std::setprecision(10) //10 is the max precision of double.
                  << double(s)
                  << std::setprecision(6); //6 is the default precision.
    }

    template<ScalarOption option, bool errorTrack>
    inline Scalar<option, errorTrack> operator+(const Scalar<option, errorTrack>& s) {
        return Scalar<option, errorTrack>(s);
    }

    template<ScalarOption option, bool errorTrack1, bool errorTrack2>
    inline void operator+=(Scalar<option, errorTrack1>& s1
            , const Scalar<option, errorTrack2>& s2) { s1 = s1 + s2; }

    template<ScalarOption option, bool errorTrack1, bool errorTrack2>
    inline void operator-=(Scalar<option, errorTrack1>& s1
            , const Scalar<option, errorTrack2>& s2) { s1 = s1 - s2; }

    template<ScalarOption option, bool errorTrack1, bool errorTrack2>
    inline void operator*=(Scalar<option, errorTrack1>& s1
            , const Scalar<option, errorTrack2>& s2) { s1 = s1 * s2; }

    template<ScalarOption option, bool errorTrack1, bool errorTrack2>
    inline void operator/=(Scalar<option, errorTrack1>& s1
            , const Scalar<option, errorTrack2>& s2) { s1 = s1 / s2; }

    template<ScalarOption option, bool errorTrack>
    inline void operator^=(Scalar<option, errorTrack>& s1
            , const Scalar<option, errorTrack>& s2) { s1 = s1 ^ s2; }

    template<ScalarOption option, bool errorTrack>
    inline void operator<<=(Scalar<option, errorTrack>& s
            , int bits) { s = s << bits; }

    template<ScalarOption option, bool errorTrack>
    inline void operator>>=(Scalar<option, errorTrack>& s
            , int bits) { s = s >> bits; }

    template<ScalarOption option>
    inline bool operator>=(const Internal::AbstractScalar<option>& s1, const Internal::AbstractScalar<option>& s2) {
        return !(s1 < s2);
    }

    template<ScalarOption option>
    inline bool operator<=(const Internal::AbstractScalar<option>& s1, const Internal::AbstractScalar<option>& s2) {
        return !(s1 > s2);
    }

    template<ScalarOption option>
    inline bool operator!= (const Internal::AbstractScalar<option>& s1, const Internal::AbstractScalar<option>& s2) {
        return !(s1 == s2);
    }
    /*!
     * The following two functions handle swap(s1, s2). Use swap of Scalars whose errorTrack is false by default.
     */
    template<ScalarOption option, bool errorTrack>
    inline void swap(Scalar<option, false>& s1, Scalar<option, errorTrack>& s2) noexcept {
        s1.swap(s2);
    }

    template<ScalarOption option, bool errorTrack>
    inline void swap(Scalar<option, true>& s1, Scalar<option, errorTrack>& s2) noexcept {
        s2.swap(s1);
    }
    ///////////////////////////////////////////MultiPrecision/////////////////////////////////////////
    template<bool errorTrack>
    inline Scalar<MultiPrecision, errorTrack>& operator++(Scalar<MultiPrecision, errorTrack>& s) {
        s += BasicConst::getInstance()._1;
        return s;
    }
    
    template<bool errorTrack>
    inline Scalar<MultiPrecision, errorTrack>& operator--(Scalar<MultiPrecision, errorTrack>& s) {
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

    inline bool absCompare(const Internal::AbstractScalar<Float>& s1, const Internal::AbstractScalar<Float>& s2) {
        return fabsf(s1.getTrivial()) >= fabsf(s2.getTrivial());
    }

    inline bool operator> (const Internal::AbstractScalar<Float>& s1, const Internal::AbstractScalar<Float>& s2) {
        return s1.getTrivial() > s2.getTrivial();
    }

    inline bool operator< (const Internal::AbstractScalar<Float>& s1, const Internal::AbstractScalar<Float>& s2) {
        return s1.getTrivial() < s2.getTrivial();
    }

    inline bool operator== (const Internal::AbstractScalar<Float>& s1, const Internal::AbstractScalar<Float>& s2) {
        return s1.getTrivial() == s2.getTrivial();
    }
    ////////////////////////////////////////Float-WithoutError///////////////////////////////////////////
    inline Scalar<Float, false>::Scalar(const Scalar<Float, true>& s) : Base(s) {}

    inline Scalar<Float, false>::Scalar(const Scalar<Double, false>& s) : Base(float(s)) {}

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
        Base::toInteger();
    }
    /////////////////////////////////////////Float-WithError////////////////////////////////////////////////
    inline void Scalar<Float, true>::toInteger() {
        Base::toInteger();
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

    inline bool absCompare(const Internal::AbstractScalar<Double>& s1, const Internal::AbstractScalar<Double>& s2) {
        return fabs(s1.getTrivial()) >= fabs(s2.getTrivial());
    }

    inline bool operator> (const Internal::AbstractScalar<Double>& s1, const Internal::AbstractScalar<Double>& s2) {
        return s1.getTrivial() > s2.getTrivial();
    }

    inline bool operator< (const Internal::AbstractScalar<Double>& s1, const Internal::AbstractScalar<Double>& s2) {
        return s1.getTrivial() < s2.getTrivial();
    }

    inline bool operator== (const Internal::AbstractScalar<Double>& s1, const Internal::AbstractScalar<Double>& s2) {
        return s1.getTrivial() == s2.getTrivial();
    }
    ////////////////////////////////////////Double-WithoutError///////////////////////////////////////////
    inline Scalar<Double, false>::Scalar(const Scalar<Float, false>& s) : Base(double(s)) {}

    inline Scalar<Double, false>::Scalar(const Scalar<Double, true>& s) : Base(s) {}

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
        Base::toInteger();
    }
    /////////////////////////////////////////Double-WithError///////////////////////////////////////////
    inline void Scalar<Double, true>::toInteger() {
        Base::toInteger();
        a = 0;
    }
}

#include "ScalarArithmetic.h"
#include "ElementaryFunction.h"
#include "Operation/Pow.h"

namespace Physica::Core {
    template<ScalarOption option>
    inline Scalar<option, false> operator^(
            const Scalar<option, false>& s1, const Scalar<option, false>& s2) {
        return Scalar<option, false>(std::pow(s1.getTrivial(), s2.getTrivial()));
    }

    template<ScalarOption option>
    inline Scalar<option, true> operator^(
            const Scalar<option, true>& s1, const Scalar<option, false>& s2) {
        const auto result = std::pow(s1.getTrivial(), s2.getTrivial());
        return Scalar<option, true>(result, std::pow(s1.getTrivial(), s2.getTrivial() + s2.getA()) - result);
    }

    template<ScalarOption option>
    inline Scalar<option, true> operator^(
            const Scalar<option, false>& s1, const Scalar<option, true>& s2) {
        const auto result = std::pow(s1.getTrivial(), s2.getTrivial());
        return Scalar<option, true>(result
                , std::pow(s1.getTrivial() + (s1.getTrivial() > 1 ? s1.getA() : -s1.getA()), s2.getTrivial()) - result);
    }

    template<ScalarOption option>
    inline Scalar<option, true> operator^(
            const Scalar<option, true>& s1, const Scalar<option, true>& s2) {
        const auto result = std::pow(s1.getTrivial(), s2.getTrivial());
        return Scalar<option, true>(result
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
