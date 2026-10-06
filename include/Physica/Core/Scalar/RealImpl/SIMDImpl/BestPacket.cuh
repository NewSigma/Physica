/*
 * Copyright 2026 Weibo He.
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

#include "BestPacket.h"
#include "Half2.h"

namespace Physica {
    template<Scalar T, size_t Length>
    class device_obj<BestPacket<T, Length>> {
        static_assert(!T::isComplex(), "[Error]: This specialization does not handle complex");
        static_assert(!T::isForwardDiff(), "[Error]: This specialization does not handle forward diff");
        static_assert(T::Prec != FloatMP, "[Error]: FloatMP is not supported on device");

        constexpr static bool isDynamic = Length == Dynamic;
        constexpr static int PackSize = []() consteval static noexcept -> int {
            if constexpr (T::Prec == Float16)
                return 2;
            else if constexpr (T::Prec == Float32)
                return 4;
            else {
                static_assert(T::Prec == Float64);
                return 2;
            }
        }();
        constexpr static bool UsePack = isDynamic || Length >= PackSize;
    public:
        constexpr static int Size = UsePack ? PackSize : 1;
        using Type = std::conditional<Size == 1, T, device_obj<SIMD<T, Size>>>::type;
    };
}
