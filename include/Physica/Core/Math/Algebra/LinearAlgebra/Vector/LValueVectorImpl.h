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
        struct InnerDotImpl {
            using ResultType = typename Internal::BinaryScalarOpReturnType<typename T1::ScalarType, typename T2::ScalarType>::Type;
            static ResultType run(const RValueVector<T1>& v1, const RValueVector<T2>& v2) {
                auto result = ResultType::Zero();
                for(size_t i = 0; i < v1.getLength(); ++i)
                    result += v1.calc(i) * v2.calc(i);
                return result;
            }
        };

        template<class T1, class T2>
        class InnerDotImpl<T1, T2, true> {
            constexpr static size_t size1 = T1::SizeAtCompile;
            constexpr static size_t size2 = T2::SizeAtCompile;
            constexpr static size_t SizeAtCompile = size1 > size2 ? size1 : size2;
        public:
            using ResultType = typename T1::ScalarType;
            using PacketType = typename Internal::BestPacket<ResultType, SizeAtCompile>::Type;

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
    LValueVector<Derived>& LValueVector<Derived>::operator=(const LValueVector& v) {
        return Base::getDerived() = v.getDerived();
    }

    template<class Derived>
    LValueVector<Derived>& LValueVector<Derived>::operator=(LValueVector&& v) noexcept {
        return Base::getDerived() = std::move(v.getDerived());
    }

    template<class Derived>
    template<class OtherVector>
    Derived& LValueVector<Derived>::operator=(const RValueVector<OtherVector>& v) {
        Base::getDerived().resize(v.getLength());
        v.getDerived().assignTo(*this);
        return Base::getDerived();
    }

    template<class Derived>
    template<class AnyScalar>
    Derived& LValueVector<Derived>::operator=(const ScalarBase<AnyScalar>& s) {
        for (size_t i = 0; i < Base::getLength(); ++i)
            (*this)[i] = s.getDerived();
        return Base::getDerived();
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

    template<class Derived>
    bool LValueVector<Derived>::isZero() const {
        for(size_t i = 0; i < Base::getLength(); ++i)
            if (!(*this)[i].isZero())
                return false;
        return true;
    }

    template<class VectorType1, class VectorType2>
    typename Internal::BinaryScalarOpReturnType<typename VectorType1::ScalarType, typename VectorType2::ScalarType>::Type
    operator*(const RValueVector<VectorType1>& v1, const RValueVector<VectorType2>& v2) {
        assert(v1.getLength() == v2.getLength());
        return Internal::InnerDotImpl<VectorType1, VectorType2, false>::run(v1, v2);
    }

    template<class Derived, class OtherDerived>
    void operator+=(LValueVector<Derived>& v1, const RValueVector<OtherDerived>& v2) {
        for (size_t i = 0; i < v1.getLength(); ++i)
            v1[i] = v1[i] + v2.calc(i);
    }

    template<class Derived, class OtherDerived>
    void operator-=(LValueVector<Derived>& v1, const RValueVector<OtherDerived>& v2) {
        for (size_t i = 0; i < v1.getLength(); ++i)
            v1[i] = v1[i] - v2.calc(i);
    }

    template<class VectorType, class ScalarType>
    inline void operator+=(LValueVector<VectorType>& v, const ScalarBase<ScalarType>& n) { v = v + n; }

    template<class VectorType, class ScalarType>
    inline void operator-=(LValueVector<VectorType>& v, const ScalarBase<ScalarType>& n) { v = v - n; }

    template<class VectorType, class ScalarType>
    inline void operator*=(LValueVector<VectorType>& v, const ScalarBase<ScalarType>& n) { v = v * n; }

    template<class VectorType, class ScalarType>
    inline void operator/=(LValueVector<VectorType>& v, const ScalarBase<ScalarType>& n) { v = v / n; }
}
