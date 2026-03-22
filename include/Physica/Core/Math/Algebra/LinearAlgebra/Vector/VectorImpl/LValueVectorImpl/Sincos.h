/*
 * Copyright 2020-2026 Weibo He.
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

#include "../LValueVector.h"

namespace Physica {
    template<Vector V1, Vector V2>
    auto sincos(const V1& x, V2& s, V2& c) {
        assert(x.getLength() == s.getLength() && x.getLength() == c.getLength());
        constexpr size_t SizeAtCompile = std::max(V1::SizeAtCompile, V2::SizeAtCompile);
        using ScalarType1 = V1::ScalarType;
        using ScalarType2 = V2::ScalarType;
        using ScalarType = Internal::BinaryScalarOpRtnTy<ScalarType1, ScalarType2>::Type;
        using PacketType = BestPacket<ScalarType, SizeAtCompile>::Type;
        if constexpr (ReverseDiff<ScalarType>) {
            size_t i = 0;
            auto result = co_for([&]{ return i < x.getLength(); }, [&]{ ++i; }, [&]{
                return sincos(x.calc(i), s[i], c[i]);
            });
            return result;
        }
        else if constexpr (PacketType::size() == 1) {
            for (size_t i = 0; i < x.getLength(); ++i)
                sincos(x.calc(i), s[i], c[i]);
        }
        else {
            constexpr int Size = PacketType::size();
            const size_t length = x.getLength();
            size_t i = 0;
            const size_t to = length / Size * Size;
            PacketType s_buffer, c_buffer;
            for (; i < to; i += Size) {
                sincos(x.template packet<Size>(i), s_buffer, c_buffer);
                s.writePacket(s_buffer, i);
                c.writePacket(c_buffer, i);
            }

            for (; i < length; ++i)
                sincos(x.calc(i), s[i], c[i]);
        }
    }
}
