/*
 * Copyright 2024 Weibo He.
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
    /**
     * \class BestPacket finds the best packet for a linear storage
     */
    template<class ScalarType, size_t Length>
    class BestPacket {
        static_assert(is_scalar<ScalarType>::value, "[Error]: Invalid ScalarType");
        static_assert(!ScalarType::isComplex, "[Error]: This specialization does not handle complex");
        static_assert(!ScalarType::isForwardDiff, "[Error]: This specialization does not handle forward diff");
        constexpr static bool isSinglePrec = ScalarType::Option == Float;
        constexpr static bool isDynamic = Length == Dynamic;
        constexpr static size_t size128 = isSinglePrec ? 4 : 2;
        constexpr static size_t size256 = isSinglePrec ? 8 : 4;
        constexpr static size_t size512 = isSinglePrec ? 16 : 8;
        constexpr static bool support128 = INSTRSET >= 2;
        constexpr static bool support256 = INSTRSET >= 7 && support128;
        constexpr static bool support512 = INSTRSET >= 9 && support256;
        constexpr static bool use128 = support128 && Length >= size128 && size128 != 0;
        constexpr static bool use256 = support256 && Length >= size256 && size256 != 0;
        constexpr static bool use512 = support512 && Length >= size512 && size512 != 0;
        constexpr static size_t Size1 = use128 ? size128 : 1;
        constexpr static size_t Size2 = use256 ? size256 : Size1;
        constexpr static size_t Size3 = use512 ? size512 : Size2;
        constexpr static size_t BiggestSize1 = support128 ? size128 : 1;
        constexpr static size_t BiggestSize2 = support256 ? size256 : BiggestSize1;
        constexpr static size_t BiggestSize = support512 ? size512 : BiggestSize2;
    public:
        constexpr static size_t Size = isDynamic ? BiggestSize : Size3;
        using Type = typename std::conditional<Size == 1, ScalarType, SIMD<ScalarType, Size>>::type;
    };

    template<class ScalarType, size_t Length>
    class device_obj<BestPacket<ScalarType, Length>> {
    public:
        constexpr static size_t Size = 1;
        using Type = ScalarType;
    };
}
