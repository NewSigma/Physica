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
#ifndef PHYSICA_BISECTIONMETHOD_H
#define PHYSICA_BISECTIONMETHOD_H
/*!
 * This file is part of implementations of \Scalar.
 * Do not include this header file, include Scalar.h instead.
 */
namespace Physica::Core {
    template<ScalarOption option = MultiPrecision, bool errorTrack = true>
    Scalar<option, errorTrack> bisectionMethod(
            Scalar<option, false> func(const Scalar<option, false>&)
            , const Scalar<option, false>& n
            , const Scalar<option, false>& x1, const Scalar<option, false>& x2) {
        Scalar<option, false> y1 = func(x1);
        Scalar<option, false> y2 = func(x2);
        return bisectionMethod<errorTrack>(func, n, x1, x2, y1, y2);
    }

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> bisectionMethod(
            Scalar<option, false> func(const Scalar<option, false>&)
            , const Scalar<option, false>& n
            , const Scalar<option, false>& x1, const Scalar<option, false>& x2
            , const Scalar<option, false>& y1, const Scalar<option, false>& y2) {
        if(n == y1)
            return Scalar<option, errorTrack>(x1);
        if(n == y2)
            return Scalar<option, errorTrack>(x2);

        if((Scalar<option, false>(n - y1).getLength() ^ Scalar<option, false>(n - y2).getLength()) >= 0) // NOLINT(hicpp-signed-bitwise)
            qFatal("Root is nonexistent.");

        Scalar<option, false> result = (x1 + x2) >> 1;
        Scalar<option, false> y_result(static_cast<SignedMPUnit>(1));

        Scalar<option, false> error = Scalar<option, false>(x1 - x2).toAbs() >> 1;
        Scalar<option, false> x_left(x1);
        Scalar<option, false> x_right(x2);
        Scalar<option, false> y_left(y1);

        bool delta_left_sign = Scalar<option, false>(n - y_left).getLength() > 0;
        bool delta_right_sign;
        do {
            y_result = func(result);
            delta_right_sign = Scalar<option, false>(n - y_result).getLength() > 0;

            if(delta_left_sign == delta_right_sign) {
                x_left = result;
                y_left = y_result;
                delta_left_sign = Scalar<option, false>(n - y_left).getLength() > 0;
            }
            else
                x_right = result;
            result = (x_left + x_right) >> 1;
            error >>= 1;
        } while(result.getPower() - error.getPower() < GlobalPrecision);
        return Scalar<option, errorTrack>(result).toUnitA();
    }
}

#endif