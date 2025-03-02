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
    template<Vector T1, Vector T2>
    class InnerDot {
        using This = InnerDot<T1, T2>;
        using Helper = Internal::EnableSIMD<T1, T2>;
        constexpr static size_t SizeAtCompile = Helper::SizeAtCompile;
        constexpr static bool isFastPacket1 = Traits<T1>::FastPacket;
        constexpr static bool isFastPacket2 = Traits<T2>::FastPacket;
        constexpr static bool enableSIMD = isFastPacket1 && isFastPacket2 && Helper::value;
    public:
        using ScalarType = Helper::ResultType;
        using PacketType = Helper::PacketType;
        constexpr static bool isForwardDiff = ScalarType::isForwardDiff;
        constexpr static bool isReverseDiff = ScalarType::isReverseDiff;
    private:
        const T1& v1;
        const T2& v2;
    public:
        InnerDot(const T1& v1_, const T2& v2_);
        InnerDot(const This&) = delete;
        InnerDot(This&&) noexcept = delete;
        ~InnerDot() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        auto calc() const;
        CoDiff<ScalarType> calc_base() const;
        ScalarType calc_mkl() const;
    };

    template<Vector T1, Vector T2>
    InnerDot<T1, T2>::InnerDot(const T1& v1_, const T2& v2_) : v1(v1_), v2(v2_) {
        assert(v1.getLength() == v2.getLength());
    }

    template<Vector T1, Vector T2>
    auto InnerDot<T1, T2>::calc() const {
        if constexpr (isForwardDiff) {
            if constexpr (!T1::isForwardDiff)
                return ScalarType(v2.values() * v1, v2.grads() * v1);
            else if constexpr (!T2::isForwardDiff)
                return ScalarType(v1.values() * v2, v1.grads() * v2);
            else {
                constexpr static int Order = ScalarType::Order - 1;
                return ScalarType(v1.values() * v2.values(), toDiffMaskVector<T2, Order>(v2) * v1.grads() + toDiffMaskVector<T1, Order>(v1) * v2.grads());
            }
        }
        else if constexpr (Internal::EnableMKL<T1, T2>::value)
            return calc_mkl();
        else
            return calc_base();
    }

    template<Vector T1, Vector T2>
    CoDiff<typename InnerDot<T1, T2>::ScalarType> InnerDot<T1, T2>::calc_base() const {
        if constexpr (isReverseDiff) {
            auto result = co_yield v1.values() * v2.values();
            if constexpr (ReverseDiff<T1>)
                v1.reverse(v2.values() * result.grad());
            if constexpr (ReverseDiff<T2>)
                v2.reverse(v1.values() * result.grad());
        }
        else if constexpr (enableSIMD) {
            const size_t length = v1.getLength();
            size_t i = 0;
            const size_t to = length / PacketType::size() * PacketType::size();
            PacketType buffer(0);
            for (; i < to; i += PacketType::size()) {
                PacketType p1 = v1.template packet<PacketType>(i);
                PacketType p2 = v2.template packet<PacketType>(i);
                buffer = mul_add(p1, p2, buffer);
            }
            if (to != length) {
                const size_t count = length - i;
                PacketType p1 = v1.template packetPartial<PacketType>(i, count);
                PacketType p2 = v2.template packetPartial<PacketType>(i, count);
                buffer = mul_add(p1, p2, buffer);
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

    template<Vector T1, Vector T2>
    [[nodiscard]] inline auto operator*(const T1& v1, const T2& v2) noexcept requires(!CUDA<T1> && !CUDA<T2>) {
        return InnerDot<T1, T2>(v1, v2).calc();
    }
}

#ifdef PHYSICA_MKL
    #include "InnerDot_MKL.h"
#endif
