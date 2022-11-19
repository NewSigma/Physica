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
        struct AddAssignImpl {
            static void run(LValueVector<T1>& v1, const RValueVector<T2>& v2) {
                using ScalarType = typename T1::ScalarType;
                for(size_t i = 0; i < v1.getLength(); ++i)
                    v1.getDerived()[i] += ScalarType(v2.calc(i));
            }
        };

        template<class T1, class T2>
        class AddAssignImpl<T1, T2, true> {
            static_assert(std::is_same<typename T1::ScalarType, typename T2::ScalarType>::value, "[Error]: SIMD on different scalars is not supported");
            static_assert(!T1::isComplex && !T2::isComplex, "[Error]: SIMD on complex scalars are not supported");
            constexpr static size_t size1 = T1::SizeAtCompile;
            constexpr static size_t size2 = T2::SizeAtCompile;
            constexpr static size_t SizeAtCompile = size1 > size2 ? size1 : size2;
            using PacketType = typename Internal::BestPacket<typename T1::ScalarType, SizeAtCompile>::Type;
        public:
            static void run(LValueVector<T1>& v1, const RValueVector<T2>& v2) {
                const size_t length = v1.getLength();
                size_t i = 0;
                const size_t to = length >= static_cast<size_t>(PacketType::size()) ? (length - PacketType::size()) : 0;
                for (; i < to; i += PacketType::size()) {
                    const PacketType sum = v1.getDerived().template packet<PacketType>(i) + v2.getDerived().template packet<PacketType>(i);
                    v1.getDerived().writePacket(i, sum);
                }
                const PacketType sum = v1.getDerived().template packetPartial<PacketType>(i) + v2.getDerived().template packetPartial<PacketType>(i);
                v1.getDerived().writePacketPartial(i, sum);
            }
        };
    }

    template<class Derived>
    ContinuousVector<Derived>& ContinuousVector<Derived>::operator=(const ContinuousVector<Derived>& v) {
        Base::operator=(v);
        return *this;
    }

    template<class Derived>
    ContinuousVector<Derived>& ContinuousVector<Derived>::operator=(ContinuousVector<Derived>&& v) noexcept {
        Base::operator=(std::forward<Base>(v));
        return *this;
    }

    template<class Derived>
    template<class PacketType>
    PacketType ContinuousVector<Derived>::packet(size_t index) const {
        assert(index + PacketType::size() <= Base::getLength());
        if constexpr (std::is_same_v<TrivialType, typename Internal::Traits<PacketType>::ScalarType>) {
            PacketType packet{};
            packet.load(data_ptr(index));
            return packet;
        }
        else
            return Base::template packet<PacketType>(index);
    }

    template<class Derived>
    template<class PacketType>
    PacketType ContinuousVector<Derived>::packetPartial(size_t index) const {
        if constexpr (std::is_same_v<TrivialType, typename Internal::Traits<PacketType>::ScalarType>) {
            PacketType packet{};
            const size_t count = Base::getLength() - index;
            packet.load_partial(count, data_ptr(index));
            return packet;
        }
        else
            return Base::template packetPartial<PacketType>(index);
    }

    template<class Derived>
    template<class PacketType>
    void ContinuousVector<Derived>::writePacket(size_t index, const PacketType packet) {
        if constexpr (std::is_same_v<TrivialType, typename Internal::Traits<PacketType>::ScalarType>)
            packet.store(data_ptr(index));
        else
            Base::template writePacket(index, packet);
    }

    template<class Derived>
    template<class PacketType>
    void ContinuousVector<Derived>::writePacketPartial(size_t index, const PacketType packet) {
        if constexpr (std::is_same_v<TrivialType, typename Internal::Traits<PacketType>::ScalarType>)
            packet.store_partial(Base::getLength() - index, data_ptr(index));
        else
            Base::template writePacketPartial(index, packet);
    }

    template<class Derived>
    std::ostream& operator<<(std::ostream& os, const ContinuousVector<Derived>& v) {
        using ScalarType = typename Derived::ScalarType;
        os.write(reinterpret_cast<const char*>(v.data_ptr(0)), v.getLength() * sizeof(ScalarType));
        return os;
    }

    template<class Derived>
    std::istream& operator>>(std::istream& is, ContinuousVector<Derived>& v) {
        using ScalarType = typename Derived::ScalarType;
        is.read(reinterpret_cast<char*>(v.data_ptr(0)), v.getLength() * sizeof(ScalarType));
        return is;
    }

    template<class Derived, class OtherDerived>
    void operator+=(ContinuousVector<Derived>& v1, const RValueVector<OtherDerived>& v2) {
        assert(v1.getLength() == v2.getLength());
        Internal::AddAssignImpl<Derived, OtherDerived, Internal::EnableSIMD<Derived, OtherDerived>::value>::run(v1, v2);
    }
}
