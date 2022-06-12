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

#include "FormatedVector.h"

namespace Physica::Core {
    namespace Internal {
        template<class T1, class T2, bool enableSIMD>
        struct AssignImpl {
            static void run(const RValueVector<T1>& v1, LValueVector<T2>& v2) {
                for (size_t i = 0; i < v1.getLength(); ++i)
                    v2[i] = v1.calc(i);
            }
        };

        template<class T1, class T2>
        class AssignImpl<T1, T2, true> {
            constexpr static size_t size1 = T1::SizeAtCompile;
            constexpr static size_t size2 = T2::SizeAtCompile;
            constexpr static size_t SizeAtCompile = size1 > size2 ? size1 : size2;
        public:
            using ScalarType = typename T1::ScalarType;
            using PacketType = typename Internal::BestPacket<ScalarType, SizeAtCompile>::Type;

            static void run(const RValueVector<T1>& v1, LValueVector<T2>& v2) {
                const size_t length = v1.getLength();
                if (length != 0) {
                    size_t i = 0;
                    const size_t to = length >= static_cast<size_t>(PacketType::size()) ? (length - PacketType::size()) : 0;
                    for (; i < to; i += PacketType::size())
                        v2.getDerived().writePacket(i, v1.getDerived().template packet<PacketType>(i));
                    v2.getDerived().writePacketPartial(i, v1.getDerived().template packetPartial<PacketType>(i));
                }
            }
        };
    }

    template<class Derived>
    template<class OtherDerived>
    void RValueVector<Derived>::assignTo(LValueVector<OtherDerived>& v) const {
        assert(v.getLength() == getLength());
        Internal::AssignImpl<Derived, OtherDerived, Internal::EnableSIMD<Derived, OtherDerived>::value>::run(*this, v);
    }

    template<class Derived>
    FormatedVector<Derived> RValueVector<Derived>::format() const {
        return FormatedVector<Derived>(*this);
    }

    template<class Derived>
    template<class PacketType>
    PacketType RValueVector<Derived>::packet(size_t index) const {
        PacketType packet{};
        for (int i = 0; i < PacketType::size(); ++i, ++index)
            packet.insert(i, calc(index).getTrivial());
        return packet;
    }

    template<class Derived>
    template<class PacketType>
    PacketType RValueVector<Derived>::packetPartial(size_t index) const {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
        PacketType packet{};
        for (int i = 0; index < getLength(); ++i, ++index)
            packet.insert(i, calc(index).getTrivial());
        return packet;
    #pragma GCC diagnostic pop
    }

    template<class Derived>
    typename RValueVector<Derived>::RealType RValueVector<Derived>::norm() const {
        return sqrt(squaredNorm());
    }

    template<class Derived>
    typename RValueVector<Derived>::RealType RValueVector<Derived>::squaredNorm() const {
        auto result = RealType::Zero();
        for(size_t i = 0; i < getLength(); ++i)
            result += calc(i).squaredNorm();
        return result;
    }

    template<class VectorType>
    typename RValueVector<VectorType>::ScalarType RValueVector<VectorType>::max() const {
        assert(getLength() != 0);
        ScalarType result = calc(0);
        for(size_t i = 1; i < getLength(); ++i) {
            ScalarType temp = calc(i);
            if (result < temp)
                result = temp;
        }
        return result;
    }

    template<class VectorType>
    typename RValueVector<VectorType>::ScalarType RValueVector<VectorType>::min() const {
        assert(getLength() != 0);
        ScalarType result = calc(0);
        for(size_t i = 1; i < getLength(); ++i) {
            ScalarType temp = calc(i);
            if (result > temp)
                result = temp;
        }
        return result;
    }

    template<class VectorType>
    typename RValueVector<VectorType>::ScalarType RValueVector<VectorType>::sum() const {
        assert(getLength() != 0);
        ScalarType result = ScalarType::Zero();
        for(size_t i = 0; i < getLength(); ++i)
            result += calc(i);
        return result;
    }

    template<class Derived>
    template<class OtherDerived>
    inline CrossProduct<Derived, OtherDerived>
    RValueVector<Derived>::crossProduct(const RValueVector<OtherDerived>& v) const noexcept {
        return CrossProduct(*this, v);
    }

    template<class Derived>
    template<class OtherDerived>
    typename RValueVector<Derived>::ScalarType
    RValueVector<Derived>::angleTo(const RValueVector<OtherDerived>& v) const noexcept {
        return arccos(Base::getDerived() * v.getDerived() / (norm() * v.norm()));
    }

    template<class VectorType>
    std::ostream& operator<<(std::ostream& os, const RValueVector<VectorType>& v) {
        size_t length = v.getLength();
        for (size_t i = 0; i < length; ++i)
            os << v.calc(i);
        return os;
    }
}
