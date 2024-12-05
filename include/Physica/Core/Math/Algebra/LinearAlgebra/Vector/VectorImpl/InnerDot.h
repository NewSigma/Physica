/*
 * Copyright 2023-2024 Weibo He.
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

#include "RValueVector.h"

namespace Physica::Core {
    template<Vector T1, Vector T2>
    class InnerDot {
        using This = InnerDot<T1, T2>;
        using Helper = Internal::EnableSIMD<T1, T2>;
        using ResultType = Helper::ResultType;
        using PacketType = Helper::PacketType;
        constexpr static bool isFastPacket1 = Traits<T1>::FastPacket;
        constexpr static bool isFastPacket2 = Traits<T2>::FastPacket;
        constexpr static bool enableSIMD = isFastPacket1 && isFastPacket2 && Helper::value;
        constexpr static bool isForwardDiff = T1::isForwardDiff || T2::isForwardDiff;
        constexpr static bool isReverseDiff = T1::isReverseDiff || T2::isReverseDiff;
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
        ResultType calc() const;
        ResultType calc_base() const;
        ResultType calc_mkl() const;
    };

    template<Vector T1, Vector T2>
    InnerDot<T1, T2>::InnerDot(const T1& v1_, const T2& v2_) : v1(v1_), v2(v2_) {
        assert(v1.getLength() == v2.getLength());
    }

    template<Vector T1, Vector T2>
    InnerDot<T1, T2>::ResultType InnerDot<T1, T2>::calc() const {
        if constexpr (isForwardDiff) {
            if constexpr (!T1::isForwardDiff)
                return ResultType(toValueVector(v2) * v1, toGradVector(v2) * v1);
            else if constexpr (!T2::isForwardDiff)
                return ResultType(toValueVector(v1) * v2, toGradVector(v1) * v2);
            else {
                constexpr static int Order = ResultType::Order - 1;
                return ResultType(toValueVector(v1) * toValueVector(v2), toDiffMaskVector<T2, Order>(v2) * toGradVector(v1) + toDiffMaskVector<T1, Order>(v1) * toGradVector(v2));
            }
        }
        else if constexpr (HasMKL() && is_continuous<T1>::value && is_continuous<T2>::value && !isReverseDiff)
            return calc_mkl();
        else
            return calc_base();
    }

    template<Vector T1, Vector T2>
    InnerDot<T1, T2>::ResultType InnerDot<T1, T2>::calc_base() const {
        if constexpr (enableSIMD) {
            const size_t length = v1.getLength();
            size_t i = 0;
            const size_t to = length / PacketType::size() * PacketType::size();
            PacketType buffer(0);
            if constexpr (isReverseDiff) {
                using PlainPacket = PacketType::PlainPacket;
                const auto head1 = v1[0];
                const auto head2 = v2[0];
                for (; i < to; i += PacketType::size()) {
                    const ResultType node1(head1.value_ptr() + i, head1.grad_ptr() + i);
                    const ResultType node2(head2.value_ptr() + i, head2.grad_ptr() + i);
                    PlainPacket p1{}, p2{};
                    p1.load(node1.value_ptr());
                    p2.load(node2.value_ptr());
                    buffer = mul_add(PacketType(p1, node1), PacketType(p2, node2), buffer);
                }
                if (to != length) {
                    const size_t count = length - i;
                    const ResultType node1(head1.value_ptr() + i, head1.grad_ptr() + i);
                    const ResultType node2(head2.value_ptr() + i, head2.grad_ptr() + i);
                    PlainPacket p1{}, p2{};
                    p1.load_partial(node1.value_ptr(), count);
                    p2.load_partial(node2.value_ptr(), count);
                    buffer = mul_add(PacketType(p1, node1), PacketType(p2, node2), buffer);
                }
            }
            else {
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
            }
            return buffer.sum();
        }
        else {
            auto result = ResultType(0);
            for(size_t i = 0; i < v1.getLength(); ++i)
                result += v1.calc(i) * v2.calc(i);
            return result;
        }
    }

    template<Vector T1, Vector T2>
    [[nodiscard]] inline auto operator*(const T1& v1, const T2& v2) noexcept {
        return InnerDot<T1, T2>(v1, v2).calc();
    }
}

#ifdef PHYSICA_MKL
    #include "InnerDotImpl/InnerDot_MKL.h"
#endif
