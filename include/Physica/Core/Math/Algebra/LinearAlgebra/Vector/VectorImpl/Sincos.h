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

#include "vectorclass/vectormath_trig.h"

namespace Physica::Core {
    template<class VectorType1, class VectorType2>
    void sincos(const RValueVector<VectorType1>& x, LValueVector<VectorType2>& s, LValueVector<VectorType2>& c) {
        assert(x.getLength() == s.getLength() && x.getLength() == c.getLength());
        constexpr static size_t Size1 = VectorType1::SizeAtCompile;
        constexpr static size_t Size2 = VectorType2::SizeAtCompile;
        constexpr static size_t SizeAtCompile = Size1 > Size2 ? Size1 : Size2;
        using ScalarType1 = typename VectorType1::ScalarType;
        using ScalarType2 = typename VectorType2::ScalarType;
        using ScalarType = typename Internal::BinaryScalarOpReturnType<ScalarType1, ScalarType2>::Type;
        using PacketType = typename Internal::BestPacket<ScalarType, SizeAtCompile>::Type;
        if constexpr (PacketType::size() == 1) {
            for (size_t i = 0; i < x.getLength(); ++i)
                sincos(x.getDerived().calc(i), s[i], c[i]);
        }
        else {
            const size_t length = x.getLength();
            size_t i = 0;
            const size_t to = length / PacketType::size() * PacketType::size();
            PacketType s_buffer, c_buffer;
            for (; i < to; i += PacketType::size()) {
                sincos(x.getDerived().template packet<PacketType>(i), s_buffer, c_buffer);
                s.getDerived().writePacket(i, s_buffer);
                c.getDerived().writePacket(i, c_buffer);
            }
            if (to != length) {
                const size_t count = length - i;
                sincos(x.getDerived().template packetPartial<PacketType>(i, count), s_buffer, c_buffer);
                s.getDerived().writePacketPartial(i, count, s_buffer);
                c.getDerived().writePacketPartial(i, count, c_buffer);
            }
        }

        constexpr bool isV2Continuous = std::is_base_of<ContinuousVector<VectorType2>, VectorType2>::value;
        constexpr bool isReverseDiff = ScalarType::isReverseDiff;
        if constexpr (isReverseDiff && isV2Continuous) {
            s.getDerived().makeContinuous();
            c.getDerived().makeContinuous();
        }
    }
}
