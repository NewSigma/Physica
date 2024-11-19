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

namespace Physica::Core {
    namespace Internal {
        template<Vector T1, Vector T2>
        class InnerDotImpl {
            using This = InnerDotImpl<T1, T2>;
        public:
            using ResultType = typename Internal::BinaryScalarOpReturnType<typename T1::ScalarType, typename T2::ScalarType>::Type;

            constexpr static size_t Size1 = T1::SizeAtCompile;
            constexpr static size_t Size2 = T2::SizeAtCompile;
            constexpr static size_t SizeAtCompile = Size1 > Size2 ? Size1 : Size2;
            using PacketType = typename BestPacket<ResultType, SizeAtCompile>::Type;

            constexpr static bool isFastPacket1 = Traits<T1>::FastPacket;
            constexpr static bool isFastPacket2 = Traits<T2>::FastPacket;
            constexpr static bool isSameType = std::is_same<typename T1::ScalarType, typename T2::ScalarType>::value;
            constexpr static bool isComplex = T1::isComplex || T2::isComplex;
            constexpr static bool isBadPacket = std::is_same<typename T1::ScalarType, PacketType>::value;
            constexpr static bool enableSIMD = isFastPacket1 && isFastPacket2 && isSameType && !isComplex && !isBadPacket;
            constexpr static bool isReverseDiff = T1::isReverseDiff || T2::isReverseDiff;
        public:
            inline static ResultType run(const T1& v1, const T2& v2);
        };

        template<Vector T1, Vector T2>
        inline typename InnerDotImpl<T1, T2>::ResultType
        InnerDotImpl<T1, T2>::run(const T1& v1, const T2& v2) {
            if constexpr (enableSIMD) {
                const size_t length = v1.getLength();
                size_t i = 0;
                const size_t to = length / PacketType::size() * PacketType::size();
                PacketType buffer(0);
                if constexpr (isReverseDiff) {
                    using PlainPacket = typename PacketType::PlainPacket;
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
                    result += ResultType(v1.calc(i)) * ResultType(v2.calc(i));
                return result;
            }
        }
    }

    template<Vector T1, Vector T2>
    [[nodiscard]] inline auto operator*(const T1& v1, const T2& v2) {
        assert(v1.getLength() == v2.getLength());
        return Internal::InnerDotImpl<T1, T2>::run(v1, v2);
    }
}