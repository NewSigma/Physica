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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/FormatedVector.h"

namespace Physica::Core {
    namespace Internal {
        template<class T1, class T2, bool enableSIMD>
        struct AssignImpl {
            static void run(const RValueVector<T1>& v1, LValueVector<T2>& v2) {
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

        template<class T1, class T2, bool enableSIMD> class InnerDotImpl;

        template<class T1, class T2>
        class Traits<InnerDotImpl<T1, T2, true>> {
            constexpr static size_t size1 = T1::SizeAtCompile;
            constexpr static size_t size2 = T2::SizeAtCompile;
        public:
            constexpr static size_t SizeAtCompile = size1 > size2 ? size1 : size2;
            using ResultType = typename T1::ScalarType;
            using PacketType = typename Internal::BestPacket<ResultType, SizeAtCompile>::Type;
        };

        template<class T1, class T2, bool enableSIMD>
        struct InnerDotImpl {
            using ResultType = typename Internal::BinaryScalarOpReturnType<typename T1::ScalarType, typename T2::ScalarType>::Type;
            static ResultType run(const RValueVector<T1>& v1, const RValueVector<T2>& v2) {
                auto result = ResultType::Zero();
                for(size_t i = 0; i < v1.getLength(); ++i)
                    result += ResultType(v1.calc(i)) * ResultType(v2.calc(i));
                return result;
            }
        };

        template<class T1, class T2>
        class InnerDotImpl<T1, T2, true> {
            static_assert(std::is_same<typename T1::ScalarType, typename T2::ScalarType>::value, "[Error]: SIMD on different scalars is not supported");
            static_assert(!T1::isComplex && !T2::isComplex, "[Error]: SIMD on complex scalars are not supported");
            using This = InnerDotImpl<T1, T2, true>;
        public:
            using ResultType = typename Traits<This>::ResultType;
            using PacketType = typename Traits<This>::PacketType;

            static ResultType run(const RValueVector<T1>& v1, const RValueVector<T2>& v2) {
                const size_t length = v1.getLength();
                size_t i = 0;
                const size_t to = length >= static_cast<size_t>(PacketType::size()) ? (length - PacketType::size()) : 0;
                PacketType buffer(0);
                for (; i < to; i += PacketType::size())
                    buffer += (v1.getDerived().template packet<PacketType>(i) * v2.getDerived().template packet<PacketType>(i));
                buffer += v1.getDerived().template packetPartial<PacketType>(i) * v2.getDerived().template packetPartial<PacketType>(i);
                return ResultType(horizontal_add(buffer));
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
        return sqrt(Base::getDerived().squaredNorm());
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

    template<class VectorType1, class VectorType2>
    typename Internal::BinaryScalarOpReturnType<typename VectorType1::ScalarType, typename VectorType2::ScalarType>::Type
    operator*(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2) {
        assert(v1.getLength() == v2.getLength());
        using ScalarType1 = typename VectorType1::ScalarType;
        using PacketType = typename Internal::Traits<Internal::InnerDotImpl<VectorType1, VectorType2, true>>::PacketType;
        constexpr bool isSameType = std::is_same<ScalarType1, typename VectorType2::ScalarType>::value;
        constexpr bool isNotComplex = !VectorType1::isComplex && !VectorType2::isComplex;
        constexpr bool isBadPacket = std::is_same<ScalarType1, PacketType>::value;
        [[maybe_unused]] constexpr bool enableSIMD = isSameType && isNotComplex && !isBadPacket;
        return Internal::InnerDotImpl<VectorType1, VectorType2, enableSIMD>::run(v1, v2);
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
