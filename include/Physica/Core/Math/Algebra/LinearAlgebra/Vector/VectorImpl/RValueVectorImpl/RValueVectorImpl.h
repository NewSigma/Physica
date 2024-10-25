/*
 * Copyright 2022-2024 Weibo He.
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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/FormatedVector.h>

namespace Physica::Core {
    namespace Internal {
        template<class T1, class T2, bool enableSIMD, class Executor>
        struct AssignImpl {
            inline static void run(const RValueVector<T1>& v1, LValueVector<T2>& v2) {
                using ScalarType = typename T2::ScalarType;
                Executor::parallel_for([&](size_t i) {
                    v2[i] = ScalarType(v1.calc(i));
                }, v1.getLength(), Executor::getNumThread()).wait();
            }
        };

        template<class T1, class T2, class Executor>
        class AssignImpl<T1, T2, true, Executor> {
            constexpr static size_t size1 = T1::SizeAtCompile;
            constexpr static size_t size2 = T2::SizeAtCompile;
            constexpr static size_t SizeAtCompile = size1 > size2 ? size1 : size2;
            constexpr static bool isT2Continuous = std::is_base_of<ContinuousVector<T2>, T2>::value;
        public:
            using ScalarType = typename T1::ScalarType;
            using AnyPacket = typename BestPacket<ScalarType, SizeAtCompile>::Type;
            constexpr static size_t PacketSize = AnyPacket::size();
            constexpr static bool isReverseDiff = ScalarType::isReverseDiff;

            inline static void run(const RValueVector<T1>& v1, LValueVector<T2>& v2) { //FIXME: The function declaration is not compatible to AddAssignImpl
                if constexpr (SizeAtCompile != Dynamic) {
                    constexpr size_t to = SizeAtCompile / PacketSize * PacketSize;
                    for (size_t i = 0; i < to; i += PacketSize)
                        v2.getDerived().writePacket(i, v1.getDerived().template packet<AnyPacket>(i));
                    
                    constexpr size_t i = SizeAtCompile - SizeAtCompile % PacketSize;
                    if constexpr (i != SizeAtCompile) {
                        constexpr size_t count = SizeAtCompile - i;
                        v2.getDerived().writePacketPartial(i, count, v1.getDerived().template packetPartial<AnyPacket>(i, count));
                    }
                }
                else {
                    const size_t length = v1.getLength();
                    if (length == 0)
                        return;

                    const size_t to = length / PacketSize * PacketSize;
                    if constexpr (std::is_same<Executor, SequentialExecutor>::value) {
                        size_t i = 0;
                        for (; i < to; i += PacketSize)
                            v2.getDerived().writePacket(i, v1.getDerived().template packet<AnyPacket>(i));

                        if (to != length) {
                            const size_t count = length - i;
                            v2.getDerived().writePacketPartial(i, count, v1.getDerived().template packetPartial<AnyPacket>(i, count));
                        }
                    }
                    else {
                        const size_t numLoop = to / PacketSize;
                        auto future = Executor::parallel_for([&v1, &v2](size_t i) {
                            const size_t i1 = i * PacketSize;
                            v2.getDerived().writePacket(i1, v1.getDerived().template packet<AnyPacket>(i1));
                        }, numLoop, Executor::getNumThread());

                        if (to != length) {
                            const size_t count = length % PacketSize;
                            v2.getDerived().writePacketPartial(to, count, v1.getDerived().template packetPartial<AnyPacket>(to, count));
                        }
                        future.wait();
                    }
                }

                if constexpr (isT2Continuous && isReverseDiff)
                    v2.getDerived().makeContinuous();
            }
        };
    }

    template<class Derived>
    template<class OtherDerived, class Executor>
    inline void RValueVector<Derived>::assignTo(LValueVector<OtherDerived>& v) const {
        constexpr size_t OtherSize = Traits<OtherDerived>::SizeAtCompile;
        static_assert(SizeAtCompile == Dynamic || OtherSize == Dynamic || SizeAtCompile == OtherSize,
                "[Error]: Size mismatch between two vector");
        assert(v.getLength() == getLength() && "[Error]: Size mismatch between two vector");
        Internal::AssignImpl<Derived, OtherDerived, Internal::EnableSIMD<Derived, OtherDerived>::value, Executor>::run(*this, v);
    }

    template<class Derived>
    template<class AnyPacket>
    inline AnyPacket RValueVector<Derived>::packet(size_t index) const {
        using T = typename Traits<AnyPacket>::ScalarType;
        T buffer[AnyPacket::size()];
        for (size_t i = 0; i < AnyPacket::size(); ++i, ++index)
            buffer[i] = T(calc(index));
        if constexpr (isExpression && isReverseDiff) { //Optimize: For expression such as A + B, there is no need to create new node
            using TracerType = typename ScalarType::TracerType;
            TracerType::getInstance().reserve(AnyPacket::size());
            for (auto& elem : buffer)
                elem = elem.copy();
        }
        AnyPacket packet{};
        packet.load(buffer);
        return packet;
    }

    template<class Derived>
    template<class AnyPacket>
    inline AnyPacket RValueVector<Derived>::packetPartial(size_t index, size_t count) const {
        using T = typename Traits<AnyPacket>::ScalarType;
        T buffer[AnyPacket::size()];
        for (size_t i = 0; i < AnyPacket::size(); ++i, ++index)
            buffer[i] = i < count ? T(calc(index)) : T(0);
        if constexpr (isExpression && isReverseDiff) {
            using TracerType = typename ScalarType::TracerType;
            TracerType::getInstance().reserve(count);
            for (size_t i = 0; i < count; ++i)
                buffer[i] = buffer[i].copy();
        }
        AnyPacket packet{};
        packet.load_partial(count, buffer);
        return packet;
    }

    template<class Derived>
    inline FormatedVector<Derived> RValueVector<Derived>::format() const {
        return FormatedVector<Derived>(*this);
    }

    template<class Derived>
    inline typename RValueVector<Derived>::RealType RValueVector<Derived>::norm1() const {
        return abs(Base::getDerived()).sum();
    }

    template<class Derived>
    inline typename RValueVector<Derived>::RealType RValueVector<Derived>::norm2() const {
        return norm();
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
    typename RValueVector<Derived>::RealType RValueVector<Derived>::lnSquaredNorm() const {
        const auto& derived = Base::getDerived();
        const RealType maxabs = abs(derived).max();
        assert(maxabs > RealType(std::numeric_limits<ScalarType>::min()) && "[Error]: ln(0) is not allowed");
        const RealType factor = reciprocal(maxabs);
        return ln((derived * factor).squaredNorm()) + RealType(2) * ln(maxabs);
    }

    template<class Derived>
    inline typename RValueVector<Derived>::RealType RValueVector<Derived>::normInf() const {
        return abs(Base::getDerived()).max();
    }

    template<class Derived>
    typename RValueVector<Derived>::ScalarType RValueVector<Derived>::max() const {
        assert(getLength() != 0);
        constexpr bool EnableSIMD = Internal::EnableSIMD<Derived>::value;
        if constexpr (!EnableSIMD || isReverseDiff) {  // Optimize: Not implemented for reverse diff
            ScalarType result = calc(0);
            for(size_t i = 1; i < getLength(); ++i) {
                ScalarType temp = calc(i);
                if (result < temp)
                    result = temp;
            }
            return result;
        }
        else {
            const auto& v = Base::getDerived();
            PacketType buffer(std::numeric_limits<ScalarType>::lowest());
            ScalarType result;
            if constexpr (SizeAtCompile != Dynamic) {
                constexpr size_t to = SizeAtCompile / PacketType::size() * PacketType::size();
                for (size_t i = 0; i < to; i += PacketType::size())
                    buffer = std::max(v.template packet<PacketType>(i), buffer);
                result = buffer.max();

                constexpr size_t i = SizeAtCompile - SizeAtCompile % PacketType::size();
                if constexpr (i != SizeAtCompile) {
                    constexpr size_t count = SizeAtCompile - i;
                    for (size_t j = 0; j < count; ++j)
                        result = std::max(result, v.calc(i + j));
                }
            }
            else {
                const size_t length = getLength();
                size_t i = 0;
                const size_t to = length / PacketType::size() * PacketType::size();
                for (; i < to; i += PacketType::size())
                    buffer = std::max(v.template packet<PacketType>(i), buffer);
                result = buffer.max();

                if (to != length) {
                    const size_t count = length - i;
                    for (size_t j = 0; j < count; ++j)
                        result = std::max(result, v.calc(i + j));
                }
            }
            return result;
        }
    }

    template<class Derived>
    typename RValueVector<Derived>::ScalarType RValueVector<Derived>::min() const {
        assert(getLength() != 0);
        constexpr bool EnableSIMD = Internal::EnableSIMD<Derived>::value;
        if constexpr (!EnableSIMD || isReverseDiff) { // Optimize: Not implemented for reverse diff
            ScalarType result = calc(0);
            for(size_t i = 1; i < getLength(); ++i) {
                ScalarType temp = calc(i);
                if (result > temp)
                    result = temp;
            }
            return result;
        }
        else {
            const auto& v = Base::getDerived();
            PacketType buffer(std::numeric_limits<ScalarType>::max());
            ScalarType result;
            if constexpr (SizeAtCompile != Dynamic) {
                constexpr size_t to = SizeAtCompile / PacketType::size() * PacketType::size();
                for (size_t i = 0; i < to; i += PacketType::size())
                    buffer = std::min(v.template packet<PacketType>(i), buffer);
                result = buffer.min();

                constexpr size_t i = SizeAtCompile - SizeAtCompile % PacketType::size();
                if constexpr (i != SizeAtCompile) {
                    constexpr size_t count = SizeAtCompile - i;
                    for (size_t j = 0; j < count; ++j)
                        result = std::min(result, v.calc(i + j));
                }
            }
            else {
                const size_t length = getLength();
                size_t i = 0;
                const size_t to = length / PacketType::size() * PacketType::size();
                for (; i < to; i += PacketType::size())
                    buffer = std::min(v.template packet<PacketType>(i), buffer);
                result = buffer.min();

                if (to != length) {
                    const size_t count = length - i;
                    for (size_t j = 0; j < count; ++j)
                        result = std::min(result, v.calc(i + j));
                }
            }
            return result;
        }
    }

    template<class Derived>
    typename RValueVector<Derived>::ScalarType RValueVector<Derived>::sum() const {
        assert(getLength() != 0);
        constexpr bool EnableSIMD = Internal::EnableSIMD<Derived>::value;
        if constexpr (!EnableSIMD || isReverseDiff) { // Optimize: Not implemented for reverse diff
            auto result = ScalarType(0);
            for(size_t i = 0; i < getLength(); ++i)
                result += calc(i);
            return result;
        }
        else {
            const auto& v = Base::getDerived();
            PacketType buffer(0);
            if constexpr (SizeAtCompile != Dynamic) {
                constexpr size_t to = SizeAtCompile / PacketType::size() * PacketType::size();
                for (size_t i = 0; i < to; i += PacketType::size())
                    buffer += v.template packet<PacketType>(i);

                constexpr size_t i = SizeAtCompile - SizeAtCompile % PacketType::size();
                if constexpr (i != SizeAtCompile) {
                    constexpr size_t count = SizeAtCompile - i;
                    buffer += v.template packetPartial<PacketType>(i, count);
                }
            }
            else {
                const size_t length = getLength();
                size_t i = 0;
                const size_t to = length / PacketType::size() * PacketType::size();
                for (; i < to; i += PacketType::size())
                    buffer += v.template packet<PacketType>(i);

                if (to != length) {
                    const size_t count = length - i;
                    buffer += v.template packetPartial<PacketType>(i, count);
                }
            }
            return buffer.sum();
        }
    }

    template<class Derived>
    typename RValueVector<Derived>::ScalarType RValueVector<Derived>::prod() const {
        assert(getLength() != 0);
        auto result = calc(0);
        for(size_t i = 1; i < getLength(); ++i)
            result *= calc(i);
        return result;
    }

    template<class Derived>
    bool RValueVector<Derived>::isZeros() const {
        for(size_t i = 0; i < getLength(); ++i)
            if (!calc(i).isZero())
                return false;
        return true;
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
    template<size_t Length>
    inline RVectorBlock<Derived, Length> RValueVector<Derived>::head(size_t to) {
        return {Base::getDerived(), 0, to};
    }

    template<class Derived>
    template<size_t Length>
    inline const RVectorBlock<Derived, Length> RValueVector<Derived>::head(size_t to) const {
        return {Base::getConstCastDerived(), 0, to};
    }

    template<class Derived>
    template<size_t Length>
    inline RVectorBlock<Derived, Length> RValueVector<Derived>::tail(size_t from) {
        return {Base::getDerived(), from};
    }

    template<class Derived>
    template<size_t Length>
    inline const RVectorBlock<Derived, Length> RValueVector<Derived>::tail(size_t from) const {
        return {Base::getConstCastDerived(), from};
    }

    template<class Derived>
    template<size_t Length>
    inline RVectorBlock<Derived, Length> RValueVector<Derived>::segment(size_t from, size_t to) {
        return {Base::getDerived(), from, to};
    }

    template<class Derived>
    template<size_t Length>
    inline const RVectorBlock<Derived, Length> RValueVector<Derived>::segment(size_t from, size_t to) const {
        return {Base::getConstCastDerived(), from, to};
    }

    template<class Derived>
    inline ReverseVector<Derived> RValueVector<Derived>::reverse() {
        return ReverseVector<Derived>(Base::getDerived());
    }

    template<class Derived>
    inline const ReverseVector<Derived> RValueVector<Derived>::reverse() const {
        return ReverseVector<Derived>(Base::getDerived());
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
}
