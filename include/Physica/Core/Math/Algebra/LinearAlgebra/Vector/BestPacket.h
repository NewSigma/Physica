/*
 * Copyright 2022 WeiBo He.
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

#include "Packet.h"
#include "Physica/Utils/Container/Array/Array.h"
#include "Physica/Core/MultiPrecision/Scalar.h"

namespace Physica::Core::Internal {
    template<class ScalarType, size_t Length>
    class PacketHelper {
        constexpr static bool isSinglePrec = ScalarType::option == Float;
        constexpr static bool isDynamic = Length == Utils::Dynamic;
        using Packet128 = typename std::conditional<isSinglePrec, Vec4f, Vec2d>::type;
        using Packet256 = typename std::conditional<isSinglePrec, Vec8f, Vec4d>::type;
        using Packet512 = typename std::conditional<isSinglePrec, Vec16f, Vec8d>::type;
        constexpr static size_t size128 = isSinglePrec ? 4 : 2;
        constexpr static size_t size256 = isSinglePrec ? 8 : 4;
        constexpr static size_t size512 = isSinglePrec ? 16 : 8;
        constexpr static bool support128 = INSTRSET >= 2;
        constexpr static bool support256 = INSTRSET >= 7;
        constexpr static bool support512 = INSTRSET >= 9;
        constexpr static bool use128 = support128 && Length >= size128;
        constexpr static bool use256 = support256 && Length >= size256;
        constexpr static bool use512 = support512 && Length >= size512;

        using SupportType1 = typename std::conditional<support128, Packet128, ScalarType>::type;
        using SupportType2 = typename std::conditional<support256, Packet256, SupportType1>::type;
        using SupportType3 = typename std::conditional<support512, Packet512, SupportType2>::type;
        using UseType1 = typename std::conditional<use128, Packet128, ScalarType>::type;
        using UseType2 = typename std::conditional<use256, Packet256, UseType1>::type;
        using UseType3 = typename std::conditional<use512, Packet512, UseType2>::type;
    public:
        using BiggestPacket = SupportType3;
        using Type = typename std::conditional<isDynamic, BiggestPacket, UseType3>::type;
        constexpr static size_t BiggestSize = support512 ? size512
                                                         : (support256 ? size256
                                                                       : (support128 ? size128 : 1));
        constexpr static size_t Size = isDynamic ? BiggestSize
                                                 : (use512 ? size512
                                                           : (use256 ? size256
                                                                     : (use128 ? size128 : 1)));
    };
    /**
     * Find the best packet for a linear storage
     */
    template<class ScalarType, size_t Length>
    class BestPacket {
        static_assert((ScalarType::option == Float || ScalarType::option == Double) && !ScalarType::errorTrack, "Unsupported float type");
    public:
        using Type = typename PacketHelper<ScalarType, Length>::Type;
        constexpr static size_t Size = PacketHelper<ScalarType, Length>::Size;
    };

    template<class ScalarType>
    struct EnableSIMDHelper {
        constexpr static bool value = false;
    };

    template<ScalarOption option>
    struct EnableSIMDHelper<Scalar<option, false>> {
        constexpr static bool value = true;
    };

    template<class VectorType1, class VectorType2>
    class EnableSIMD {
        constexpr static bool same_scalar = std::is_same<typename VectorType1::ScalarType, typename VectorType2::ScalarType>::value;
    public:
        constexpr static bool value = same_scalar && EnableSIMDHelper<typename VectorType1::ScalarType>::value;
    };
}
