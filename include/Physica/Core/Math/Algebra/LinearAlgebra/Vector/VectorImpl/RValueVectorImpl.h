/*
 * Copyright 2022-2023 WeiBo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/FormatedVector.h"

namespace Physica::Core {
    namespace Internal {
        template<class T1, class T2, bool enableSIMD>
        struct AssignImpl {
            inline static void run(const RValueVector<T1>& v1, LValueVector<T2>& v2) {
                using ScalarType = typename T2::ScalarType;
                for (size_t i = 0; i < v1.getLength(); ++i)
                    v2[i] = ScalarType(v1.calc(i));
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

            inline static void run(const RValueVector<T1>& v1, LValueVector<T2>& v2) {
                if constexpr (SizeAtCompile != Dynamic) {
                    constexpr size_t to = SizeAtCompile / PacketType::size() * PacketType::size();
                    for (size_t i = 0; i < to; i += PacketType::size())
                        v2.getDerived().writePacket(i, v1.getDerived().template packet<PacketType>(i));
                    
                    constexpr size_t i = SizeAtCompile - SizeAtCompile % PacketType::size();
                    if constexpr (i != SizeAtCompile) {
                        constexpr size_t count = SizeAtCompile - i;
                        v2.getDerived().writePacketPartial(i, count, v1.getDerived().template packetPartial<PacketType>(i, count));
                    }
                }
                else {
                    const size_t length = v1.getLength();
                    if (length != 0) {
                        size_t i = 0;
                        const size_t to = length / PacketType::size() * PacketType::size();
                        for (; i < to; i += PacketType::size())
                            v2.getDerived().writePacket(i, v1.getDerived().template packet<PacketType>(i));
                        if (to != length) {
                            const size_t count = length - i;
                            v2.getDerived().writePacketPartial(i, count, v1.getDerived().template packetPartial<PacketType>(i, count));
                        }
                    }
                }
            }
        };

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
        public:
            inline static ResultType run(const RValueVector<T1>& v1, const RValueVector<T2>& v2) {
                if constexpr (enableSIMD) {
                    const size_t length = v1.getLength();
                    size_t i = 0;
                    const size_t to = length / PacketType::size() * PacketType::size();
                    PacketType buffer(0);
                    for (; i < to; i += PacketType::size())
                        buffer += (v1.getDerived().template packet<PacketType>(i) * v2.getDerived().template packet<PacketType>(i));
                    if (to != length) {
                        const size_t count = length - i;
                        buffer += v1.getDerived().template packetPartial<PacketType>(i, count) * v2.getDerived().template packetPartial<PacketType>(i, count);
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
        };
    }

    template<class Derived>
    template<class OtherDerived>
    void RValueVector<Derived>::assignTo(LValueVector<OtherDerived>& v) const {
        constexpr size_t OtherSize = Internal::Traits<OtherDerived>::SizeAtCompile;
        static_assert(SizeAtCompile == Dynamic || OtherSize == Dynamic || SizeAtCompile == OtherSize,
                "[Error]: Size mismatch between two vector");
        assert(v.getLength() == getLength());
        Internal::AssignImpl<Derived, OtherDerived, Internal::EnableSIMD<Derived, OtherDerived>::value>::run(*this, v);
    }

    template<class Derived>
    FormatedVector<Derived> RValueVector<Derived>::format() const {
        return FormatedVector<Derived>(*this);
    }

    template<class Derived>
    template<class PacketType>
    inline PacketType RValueVector<Derived>::packet(size_t index) const {
        PacketType packet{};
        for (size_t i = 0; i < PacketType::size(); ++i, ++index)
            packet.insert(i, calc(index).getTrivial());
        return packet;
    }

    template<class Derived>
    template<class PacketType>
    inline PacketType RValueVector<Derived>::packetPartial(size_t index, size_t count) const {
        auto packet = PacketType(0);
        for (size_t i = 0; i < count; ++i, ++index)
            packet.insert(i, calc(index).getTrivial());
        return packet;
    }

    template<class Derived>
    inline typename RValueVector<Derived>::RealType RValueVector<Derived>::norm() const {
        return sqrt(Base::getDerived().squaredNorm());
    }

    template<class Derived>
    inline typename RValueVector<Derived>::RealType RValueVector<Derived>::squaredNorm() const {
        auto result = RealType(0);
        for(size_t i = 0; i < getLength(); ++i)
            result += calc(i).squaredNorm();
        return result;
    }

    template<class Derived>
    typename RValueVector<Derived>::ScalarType RValueVector<Derived>::max() const {
        assert(getLength() != 0);
        ScalarType result = calc(0);
        for(size_t i = 1; i < getLength(); ++i) {
            ScalarType temp = calc(i);
            if (result < temp)
                result = temp;
        }
        return result;
    }

    template<class Derived>
    typename RValueVector<Derived>::ScalarType RValueVector<Derived>::min() const {
        assert(getLength() != 0);
        ScalarType result = calc(0);
        for(size_t i = 1; i < getLength(); ++i) {
            ScalarType temp = calc(i);
            if (result > temp)
                result = temp;
        }
        return result;
    }

    template<class Derived>
    typename RValueVector<Derived>::ScalarType RValueVector<Derived>::sum() const {
        assert(getLength() != 0);
        auto result = ScalarType(0);
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

    template<class Derived>
    inline RVectorBlock<Derived> RValueVector<Derived>::head(size_t to) {
        return RVectorBlock<Derived>(Base::getDerived(), 0, to);
    }

    template<class Derived>
    inline const RVectorBlock<Derived> RValueVector<Derived>::head(size_t to) const {
        return RVectorBlock<Derived>(Base::getConstCastDerived(), 0, to);
    }

    template<class Derived>
    inline RVectorBlock<Derived> RValueVector<Derived>::tail(size_t from) {
        return RVectorBlock<Derived>(Base::getDerived(), from);
    }

    template<class Derived>
    inline const RVectorBlock<Derived> RValueVector<Derived>::tail(size_t from) const {
        return RVectorBlock<Derived>(Base::getConstCastDerived(), from);
    }

    template<class Derived>
    inline RVectorBlock<Derived> RValueVector<Derived>::segment(size_t from, size_t to) {
        return RVectorBlock<Derived>(Base::getDerived(), from, to);
    }

    template<class Derived>
    inline const RVectorBlock<Derived> RValueVector<Derived>::segment(size_t from, size_t to) const {
        return RVectorBlock<Derived>(Base::getConstCastDerived(), from, to);
    }

    template<class Derived>
    inline ReverseVector<Derived> RValueVector<Derived>::reverse() {
        return ReverseVector<Derived>(Base::getDerived());
    }

    template<class Derived>
    inline const ReverseVector<Derived> RValueVector<Derived>::reverse() const {
        return ReverseVector<Derived>(Base::getConstDerived());
    }

    template<class VectorType1, class VectorType2>
    typename Internal::BinaryScalarOpReturnType<typename VectorType1::ScalarType, typename VectorType2::ScalarType>::Type
    operator*(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2) {
        assert(v1.getLength() == v2.getLength());
        return Internal::InnerDotImpl<VectorType1, VectorType2>::run(v1, v2);
    }

    template<class VectorType1, class VectorType2>
    bool vectorNear(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2, double precision) {
        using ScalarType = typename Internal::BinaryScalarOpReturnType<typename VectorType1::ScalarType, typename VectorType2::ScalarType>::Type;
        assert(v1.getLength() == v2.getLength());
        for (size_t i = 0; i < v1.getLength(); ++i)
            if (!scalarNear(ScalarType(v1.calc(i)), ScalarType(v2.calc(i)), precision))
                return false;
        return true;
    }

    template<class VectorType>
    std::ostream& operator<<(std::ostream& os, const RValueVector<VectorType>& v) {
        os << v.getDerived();
        return os;
    }
}
