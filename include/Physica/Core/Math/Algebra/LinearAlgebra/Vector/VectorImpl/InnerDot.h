/*
 * Copyright 2023 WeiBo He.
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
        template<class T1, class T2>
        class InnerDotImpl {
            using This = InnerDotImpl<T1, T2>;
        public:
            using ResultType = typename Internal::BinaryScalarOpReturnType<typename T1::ScalarType, typename T2::ScalarType>::Type;

            constexpr static size_t size1 = T1::SizeAtCompile;
            constexpr static size_t size2 = T2::SizeAtCompile;
            constexpr static size_t SizeAtCompile = size1 > size2 ? size1 : size2;
            using PacketType = typename Internal::BestPacket<ResultType, SizeAtCompile>::Type;

            constexpr static bool isSameType = std::is_same<typename T1::ScalarType, typename T2::ScalarType>::value;
            constexpr static bool isNotComplex = !T1::isComplex && !T2::isComplex;
            constexpr static bool isBadPacket = std::is_same<typename T1::ScalarType, PacketType>::value;
            constexpr static bool enableSIMD = isSameType && isNotComplex && !isBadPacket;

            constexpr static bool isT1Continuous = std::is_base_of<ContinuousVector<T1>, T1>::value;
            constexpr static bool isT2Continuous = std::is_base_of<ContinuousVector<T2>, T2>::value;
            constexpr static bool isReverseDiff = T1::isReverseDiff || T2::isReverseDiff;
        public:
            inline static ResultType run(const RValueVector<T1>& v1, const RValueVector<T2>& v2);
        };

        template<class T1, class T2>
        inline typename InnerDotImpl<T1, T2>::ResultType
        InnerDotImpl<T1, T2>::run(const RValueVector<T1>& v1, const RValueVector<T2>& v2) {
            if constexpr (enableSIMD) {
                const size_t length = v1.getLength();
                size_t i = 0;
                const size_t to = length / PacketType::size() * PacketType::size();
                PacketType buffer(0);
                if constexpr (isT1Continuous && isT2Continuous && isReverseDiff) {
                    using PlainPacket = typename PacketType::PlainPacket;
                    const auto head1 = v1.getDerived()[0];
                    const auto head2 = v2.getDerived()[0];
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
                        p1.load_partial(count, node1.value_ptr());
                        p2.load_partial(count, node2.value_ptr());
                        buffer = mul_add(PacketType(p1, node1), PacketType(p2, node2), buffer);
                    }
                }
                else {
                    for (; i < to; i += PacketType::size()) {
                        PacketType p1 = v1.getDerived().template packet<PacketType>(i);
                        PacketType p2 = v2.getDerived().template packet<PacketType>(i);
                        buffer = mul_add(p1, p2, buffer);
                    }
                    if (to != length) {
                        const size_t count = length - i;
                        PacketType p1 = v1.getDerived().template packetPartial<PacketType>(i, count);
                        PacketType p2 = v2.getDerived().template packetPartial<PacketType>(i, count);
                        buffer = mul_add(p1, p2, buffer);
                    }
                }
                return buffer.horizontal_add();
            }
            else {
                auto result = ResultType(0);
                for(size_t i = 0; i < v1.getLength(); ++i)
                    result += ResultType(v1.calc(i)) * ResultType(v2.calc(i));
                return result;
            }
        }
    }

    template<class VectorType1, class VectorType2>
    typename Internal::BinaryScalarOpReturnType<typename VectorType1::ScalarType, typename VectorType2::ScalarType>::Type
    operator*(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2) {
        assert(v1.getLength() == v2.getLength());
        return Internal::InnerDotImpl<VectorType1, VectorType2>::run(v1.getDerived(), v2.getDerived());
    }
}