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
#ifndef PHYSICA_NUMBERTHEORY_H
#define PHYSICA_NUMBERTHEORY_H

#include "Physica/Core/MultiPrecition/Scalar.h"

namespace Physica::Core {
    /*!
     * Return if s is a odd number.
     */
    template<ScalarType type>
    bool isOdd(const Scalar<type, false>& s) {
        return s.isInteger() && (static_cast<unsigned int>(s.getTrivial()) & 1U);
    }

    template<>
    bool isOdd(const Scalar<MultiPrecision, false>& s) {
        return s.isInteger() && (s[0] & 1U);
    }

    template<ScalarType type>
    bool isEven(const Scalar<type, false>& s) {
        return s.isInteger() && !(static_cast<unsigned int>(s.getTrivial()) & 1U);
    }

    template<>
    bool isEven(const Scalar<MultiPrecision, false>& s) {
        return s.isInteger() && !(s[0] & 1U);
    }
    /*!
     * Optimize: Make use of DP and the fact that $C_n^m = C_n^(n - m)$
     */
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> bernoulli(const Scalar<type, errorTrack>& s) {
        Q_ASSERT(s.isInteger() && !s.isNegative());
        Scalar<type, errorTrack> result = Scalar<type, errorTrack>::getOne();
        if(s.isZero())
            return result;
        Scalar<type, false> temp = Scalar<type, false>::getZero();
        for(; temp < s; ++temp)
            result -= combination(s, temp) * bernoulli(temp) / (s - temp + 1);
        return result;
    }
    /*!
     * Calculates the division which is a integer. s1 and s2 must be integer.
     */
    template<ScalarType type>
    Scalar<type, false> div(const Scalar<type, false>& s1, const Scalar<type, false>& s2) {
        Q_ASSERT(s1.isInteger() && s2.isInteger());
        return (s1 / s2).toInteger();
    }
    /*!
     * Calculate the remainder. s1 and s2 must be integer.
     */
    template<ScalarType type>
    Scalar<type, false> mod(const Scalar<type, false>& s1, const Scalar<type, false>& s2) {
        Q_ASSERT(s1.isInteger() && s2.isInteger());
        return s1 - div(s1, s2) * s2;
    }
    /*!
     * Implementation of Euclidean algorithm.
     * Different from decreasesTechnique(), it has better performance to big integers.
     */
    template<ScalarType type>
    Scalar<type, false> euclideanAlgorithm(const Scalar<type, false>& s1, const Scalar<type, false>& s2) {
        Q_ASSERT(s1.isInteger() && s2.isInteger());
        Scalar<type, false>* bigger, *smaller;
        if(s1 < s2) {
            smaller = s1;
            bigger = s2;
        }
        else {
            smaller = s2;
            bigger = s1;
        }
        Scalar<type, false> x(bigger), y(smaller);
        while(!y.isZero()) {
            x = mod(x, y);
            swap(x, y);
        }
        return x;
    }
    /*!
     * Implementation of decreases technique.
     * Different from euclideanAlgorithm(),
     * it has better performance to small integers because subtract is faster than division.
     *
     * It is recommended that s1 and s2 are less than INT_MAX(approximately) to avoid overflow.
     */
    template<ScalarType type>
    Scalar<type, false> decreasesTechnique(const Scalar<type, false>& s1, const Scalar<type, false>& s2) {
        Q_ASSERT(s1.isInteger() && s2.isInteger());
        Scalar<type, false>* bigger, *smaller;
        if(s1 < s2) {
            smaller = s1;
            bigger = s2;
        }
        else {
            smaller = s2;
            bigger = s1;
        }
        Scalar<type, false> x(bigger), y(smaller);
        int shift = 0; //May overflow here
        while(isEven(x) && isEven(y)) {
            Q_ASSERT(shift < INT_MAX); //Add a check
            x >>= 1;
            y >>= 1;
            ++shift;
        }

        while(!y.isZero()) {
            x = x - y;
            swap(x, y);
        }
        return x << shift;
    }
}

#endif