/*
 * Copyright 2020-2022 Weibo He.
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
    template<Vector T1, LVector T2>
    void sincos(const T1& x, T2& s, T2& c) {
        assert(x.getLength() == s.getLength() && x.getLength() == c.getLength());
        constexpr static size_t Size1 = T1::SizeAtCompile;
        constexpr static size_t Size2 = T2::SizeAtCompile;
        constexpr static size_t SizeAtCompile = Size1 > Size2 ? Size1 : Size2;
        using ScalarType1 = T1::ScalarType;
        using ScalarType2 = T2::ScalarType;
        using ScalarType = Internal::BinaryScalarOpRtnTy<ScalarType1, ScalarType2>::Type;
        using PacketType = BestPacket<ScalarType, SizeAtCompile>::Type;
        if constexpr (PacketType::size() == 1) {
            for (size_t i = 0; i < x.getLength(); ++i)
                sincos(x.calc(i), s[i], c[i]);
        }
        else {
            const size_t length = x.getLength();
            size_t i = 0;
            const size_t to = length / PacketType::size() * PacketType::size();
            PacketType s_buffer, c_buffer;
            for (; i < to; i += PacketType::size()) {
                sincos(x.template packet<PacketType>(i), s_buffer, c_buffer);
                s.writePacket(i, s_buffer);
                c.writePacket(i, c_buffer);
            }
            if (to != length) {
                const size_t count = length - i;
                sincos(x.template packetPartial<PacketType>(i, count), s_buffer, c_buffer);
                s.writePacketPartial(i, count, s_buffer);
                c.writePacketPartial(i, count, c_buffer);
            }
        }

        constexpr bool isContinuous = is_continuous<T2>::value;
        constexpr bool isReverseDiff = ScalarType::isReverseDiff;
        if constexpr (isContinuous && isReverseDiff) {
            s.makeContinuous();
            c.makeContinuous();
        }
    }
}
