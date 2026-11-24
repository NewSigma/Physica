/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
