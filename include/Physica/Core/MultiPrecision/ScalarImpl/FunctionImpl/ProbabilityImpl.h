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
    //!Return a real number between 0 and 1.
    template<class ScalarType>
    inline ScalarType randomScalar() {
        return ScalarType(static_cast<double>(random()) / RAND_MAX);
    }
    //!Return a real number between lowerBound and upperBound.
    template<class ScalarType>
    inline ScalarType randomScalar(const ScalarBase<ScalarType>& lowerBound, const ScalarBase<ScalarType>& upperBound) {
        using ResultType = typename Internal::RemoveErrorTrack<ScalarType>::Type;
        auto castUpper = ResultType(upperBound.getDerived());
        auto castLower = ResultType(lowerBound.getDerived());
        auto s = randomScalar<ResultType>();
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
    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> floor(const Scalar<option, errorTrack>& s) {
        return Scalar<option, errorTrack>(static_cast<size_t>(s.getTrivial()));
    }

    template<ScalarOption option, bool errorTrack>
    inline Scalar<option, errorTrack> ceil(const Scalar<option, errorTrack>& s) {
        return ++floor(s);
    }
}
