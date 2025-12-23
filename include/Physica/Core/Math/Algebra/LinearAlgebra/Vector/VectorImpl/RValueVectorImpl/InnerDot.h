/*
 * Copyright 2023-2025 Weibo He.
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

#include "../RValueVector.h"

namespace Physica {
    template<Vector V1, Vector V2>
    class InnerDot {
        using This = InnerDot<V1, V2>;
        using Helper = Internal::EnableSIMD<V1, V2>;
        constexpr static size_t SizeAtCompile = Helper::SizeAtCompile;
        constexpr static bool isFastPacket1 = Traits<V1>::FastPacket;
        constexpr static bool isFastPacket2 = Traits<V2>::FastPacket;
        constexpr static bool enableSIMD = isFastPacket1 && isFastPacket2 && Helper::value;
    public:
        using T1 = V1::ScalarType;
        using T2 = V2::ScalarType;
        using ScalarType = Helper::ResultType;
        using PacketType = Helper::PacketType;
        constexpr static bool isForwardDiff = ScalarType::isForwardDiff;
        constexpr static bool isReverseDiff = ScalarType::isReverseDiff;
    private:
        const V1& v1;
        const V2& v2;
    public:
        InnerDot(const V1& v1_, const V2& v2_);
        InnerDot(const This&) = delete;
        InnerDot(This&&) noexcept = delete;
        ~InnerDot() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        auto calc() const noexcept;
        CoDiff<ScalarType> calc_base() const noexcept;
        ScalarType calc_mkl() const noexcept;
    };

    template<Vector V1, Vector V2>
    InnerDot<V1, V2>::InnerDot(const V1& v1_, const V2& v2_) : v1(v1_), v2(v2_) {
        assert(v1.getLength() == v2.getLength());
    }

    template<Vector V1, Vector V2>
    auto InnerDot<V1, V2>::calc() const noexcept {
        if constexpr (isForwardDiff) {
            if constexpr (!V1::isForwardDiff)
                return ScalarType(v2.values() * v1, v2.grads() * v1);
            else if constexpr (!V2::isForwardDiff)
                return ScalarType(v1.values() * v2, v1.grads() * v2);
            else {
                constexpr static int Order = ScalarType::Order - 1;
                return ScalarType(v1.values() * v2.values(), toDiffMaskVector<V2, Order>(v2) * v1.grads() + toDiffMaskVector<V1, Order>(v1) * v2.grads());
            }
        }
        else if constexpr (Internal::EnableMKL<V1, V2>::value)
            return calc_mkl();
        else
            return calc_base();
    }

    template<Vector V1, Vector V2>
    auto InnerDot<V1, V2>::calc_base() const noexcept -> CoDiff<ScalarType> {
        if constexpr (isReverseDiff) {
            auto result = co_yield v1.values() * v2.values();
            if constexpr (ReverseDiff<V1>)
                v1.reverse(v2.values() * result.grad());
            if constexpr (ReverseDiff<V2>)
                v2.reverse(v1.values() * result.grad());
        }
        else if constexpr (enableSIMD) {
            const size_t length = v1.getLength();
            size_t i = 0;
            const size_t to = length / PacketType::size() * PacketType::size();
            auto buffer = PacketType::zeros();
            for (; i < to; i += PacketType::size()) {
                PacketType p1 = v1.template packet<PacketType>(i);
                PacketType p2 = v2.template packet<PacketType>(i);
                buffer = fma(p1, p2, buffer);
            }
            if (to != length) {
                const size_t count = length - i;
                PacketType p1 = v1.template packetPartial<PacketType>(i, count);
                PacketType p2 = v2.template packetPartial<PacketType>(i, count);
                buffer = fma(p1, p2, buffer);
            }
            co_return buffer.sum();
        }
        else {
            auto result = ScalarType(0);
            for(size_t i = 0; i < v1.getLength(); ++i)
                result += v1.calc(i) * v2.calc(i);
            co_return std::move(result);
        }
    }

    template<Vector V1, Vector V2>
    [[nodiscard]] auto operator*(const V1& v1, const V2& v2) noexcept requires(!CUDA<V1> && !CUDA<V2>) {
        return InnerDot<V1, V2>(v1, v2).calc();
    }
}

#ifdef PHYSICA_MKL
    #include "InnerDot_MKL.h"
#endif
