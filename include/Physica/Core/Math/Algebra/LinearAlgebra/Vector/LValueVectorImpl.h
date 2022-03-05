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

namespace Physica::Core {
    namespace Internal {
        template<class T1, class T2, bool enableSIMD>
        struct AssignImpl {
            static void run(const LValueVector<T1>& v1, LValueVector<T2>& v2) {
                for (size_t i = 0; i < v1.getLength(); ++i)
                    v2[i] = v1[i];
            }
        };

        template<class T1, class T2>
        class AssignImpl<T1, T2, true> {
            constexpr static size_t size1 = T1::SizeAtCompile;
            constexpr static size_t size2 = T2::SizeAtCompile;
            constexpr static size_t SizeAtCompile = size1 > size2 ? size1 : size2;
        public:
            using ScalarType = typename T1::ScalarType;
            using TrivialType = typename ScalarType::TrivialType;
            using PointerType = typename std::add_pointer<TrivialType>::type;
            using PacketType = typename Internal::BestPacket<ScalarType, SizeAtCompile>::Type;

            static void run(const LValueVector<T1>& v1, LValueVector<T2>& v2) {
                const size_t length = v1.getLength();
                if (length != 0) {
                    size_t i = 0;
                    const size_t to = length >= static_cast<size_t>(PacketType::size()) ? (length - PacketType::size()) : 0;
                    for (; i < to; i += PacketType::size())
                        v2.writePacket(i, v1.template packet<PacketType>(i));
                    v2.writePacketPartial(i, v1.template packetPartial<PacketType>(i));
                }
            }
        };
    }

    template<class Derived>
    template<class OtherDerived>
    void LValueVector<Derived>::assignTo(LValueVector<OtherDerived>& v) const {
        assert(v.getLength() == Base::getLength());
        Internal::AssignImpl<Derived, OtherDerived, Internal::EnableSIMD<Derived, OtherDerived>::value>::run(*this, v);
    }

    template<class Derived>
    template<class PacketType>
    PacketType LValueVector<Derived>::packet(size_t index) const {
        PacketType packet{};
        for (int i = 0; i < PacketType::size(); ++i, ++index)
            packet.insert(i, (*this)[index].getTrivial());
        return packet;
    }

    template<class Derived>
    template<class PacketType>
    PacketType LValueVector<Derived>::packetPartial(size_t index) const {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
        PacketType packet{};
        for (int i = 0; index < Base::getLength(); ++i, ++index)
            packet.insert(i, (*this)[index].getTrivial());
        return packet;
    #pragma GCC diagnostic pop
    }

    template<class Derived>
    template<class PacketType>
    void LValueVector<Derived>::writePacket(size_t index, const PacketType packet) {
        using TrivialType = typename ScalarType::TrivialType;
        TrivialType buffer[PacketType::size()];
        packet.store(buffer);
        for (int i = 0; i < PacketType::size(); ++i, ++index)
            (*this)[index] = buffer[i];
    }

    template<class Derived>
    template<class PacketType>
    void LValueVector<Derived>::writePacketPartial(size_t index, const PacketType packet) {
        using TrivialType = typename ScalarType::TrivialType;
        TrivialType buffer[PacketType::size()];
        packet.store(buffer);
        for (int i = 0; index < Base::getLength(); ++i, ++index)
            (*this)[index] = buffer[i];
    }
}
