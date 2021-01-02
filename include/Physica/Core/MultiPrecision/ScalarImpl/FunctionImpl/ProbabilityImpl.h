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
#ifndef PHYSICA_PROBABILITYIMPL_H
#define PHYSICA_PROBABILITYIMPL_H
/*!
 * This file is part of implementations of \Scalar.
 * Do not include this header file, include Scalar.h instead.
 */
namespace Physica::Core {
    //!Return a real number between 0 and 1.
    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> randomScalar() {
        return Scalar<type, errorTrack>(static_cast<double>(random()) / RAND_MAX);
    }
    //!Return a real number lowerBound and upperBound.
    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> randomScalar(
            const Scalar<type, errorTrack>& lowerBound,
            const Scalar<type, errorTrack>& upperBound) {
        auto castUpper = static_cast<Scalar<type, false>>(upperBound);
        auto castLower = static_cast<Scalar<type, false>>(lowerBound);
        auto s = randomScalar<type, errorTrack>();
        return s * (castUpper - castLower) + castLower;
    }

    template<bool errorTrack>
    Scalar<MultiPrecision, errorTrack> floor(const Scalar<MultiPrecision, errorTrack>& s) {
        if(s.isInteger())
            return Scalar(s);
        const auto size = s.getSize();
        const auto power = s.getPower();
        const auto power_1 = power + 1;
        auto length = size > power_1 ? power_1 : size;
        length = s.isNegative() ? -length : length;
        Scalar<MultiPrecision, errorTrack> result(length, power);
        for(int i = 0; i < length; ++i)
            result.setByte(i, s[i]);
        return result;
    }
    /*!
     * Specialization for float and double.
     * Fix: May overflow.
     */
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> floor(const Scalar<type, errorTrack>& s) {
        return Scalar<type, errorTrack>(static_cast<size_t>(s.getTrivial()));
    }

    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> ceil(const Scalar<type, errorTrack>& s) {
        return ++floor(s);
    }
}

#endif