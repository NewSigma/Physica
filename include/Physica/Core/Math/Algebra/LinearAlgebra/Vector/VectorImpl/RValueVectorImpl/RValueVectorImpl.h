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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/FormatedVector.h"

namespace Physica::Core {
    namespace Internal {
        template<LVector T1, Vector T2, class Executor>
        class AssignImpl {
            constexpr static size_t Size1 = T1::SizeAtCompile;
            constexpr static size_t Size2 = T2::SizeAtCompile;
            constexpr static size_t SizeAtCompile = Size1 > Size2 ? Size1 : Size2;
            static_assert(Size1 == Dynamic || Size2 == Dynamic || Size1 == Size2, "[Error]: Size mismatch between two vector");
            static_assert(T1::isComplex || !T2::isComplex, "[Error]: Cannot convert a complex to a real");
        public:
            using ScalarType = T1::ScalarType;
            using AnyPacket = BestPacket<ScalarType, SizeAtCompile>::Type;
            constexpr static size_t PacketSize = AnyPacket::size();

            inline static void run(T1& v1, const T2& v2) {
                assert(v1.getLength() == v2.getLength() && "[Error]: Size mismatch between two vector");
                if constexpr (Internal::EnableSIMD<T1, T2>::value)
                    run_simd(v1, v2);
                else
                    run_for(v1, v2);
            }
        private:
            inline static void run_for(T1& v1, const T2& v2) {
                using ScalarType = T1::ScalarType;
                Executor::parallel_for([&](size_t i) {
                    v1[i] = ScalarType(v2.calc(i));
                }, v2.getLength(), Executor::getNumThread()).wait();
            }

            inline static void run_simd(T1& v1, const T2& v2) {
                if constexpr (SizeAtCompile != Dynamic) {
                    constexpr size_t to = SizeAtCompile / PacketSize * PacketSize;
                    for (size_t i = 0; i < to; i += PacketSize)
                        v1.writePacket(i, v2.template packet<AnyPacket>(i));
                    
                    constexpr size_t i = SizeAtCompile - SizeAtCompile % PacketSize;
                    if constexpr (i != SizeAtCompile) {
                        constexpr size_t count = SizeAtCompile - i;
                        v1.writePacketPartial(i, count, v2.template packetPartial<AnyPacket>(i, count));
                    }
                }
                else {
                    const size_t length = v2.getLength();
                    if (length == 0)
                        return;

                    const size_t to = length / PacketSize * PacketSize;
                    if constexpr (std::is_same<Executor, SequentialExecutor>::value) {
                        size_t i = 0;
                        for (; i < to; i += PacketSize)
                            v1.writePacket(i, v2.template packet<AnyPacket>(i));

                        if (to != length) {
                            const size_t count = length - i;
                            v1.writePacketPartial(i, count, v2.template packetPartial<AnyPacket>(i, count));
                        }
                    }
                    else {
                        const size_t numLoop = to / PacketSize;
                        auto future = Executor::parallel_for([&v2, &v1](size_t i) {
                            const size_t i1 = i * PacketSize;
                            v1.writePacket(i1, v2.template packet<AnyPacket>(i1));
                        }, numLoop, Executor::getNumThread());

                        if (to != length) {
                            const size_t count = length % PacketSize;
                            v1.writePacketPartial(to, count, v2.template packetPartial<AnyPacket>(to, count));
                        }
                        future.wait();
                    }
                }
            }
        };
    }

    template<class Derived>
    template<LVector V, class Executor>
    inline void RValueVector<Derived>::assignTo(V& v) const {
        Internal::AssignImpl<V, Derived, Executor>::run(v, Base::getDerived());
    }

    template<class Derived>
    template<class AnyPacket>
    inline AnyPacket RValueVector<Derived>::packet(size_t index) const {
        using T = Traits<AnyPacket>::ScalarType;
        assert(index + AnyPacket::size() <= getLength() && "[Error]: Index out of range");
        if constexpr (Diffable<T>) {
            using ValuePacket = AnyPacket::ValueType;
            if constexpr (isForwardDiff) {
                using GradPacket = AnyPacket::GradType;
                const auto& x = Base::getDerived();
                auto values = toValueVector(x).template packet<ValuePacket>(index);
                auto grads = toGradVector(x).template packet<GradPacket>(index);
                return AnyPacket(std::move(values), std::move(grads));
            }
            else {
                using ValueType = T::ValueType;
                ValueType values[AnyPacket::size()];
                for (size_t i = 0; i < AnyPacket::size(); ++i, ++index)
                    values[i] = ValueType(calc(index));
                ValuePacket packet{};
                packet.load(values);
                return AnyPacket(std::move(packet));
            }
        }
        else {
            T buffer[AnyPacket::size()];
            for (size_t i = 0; i < AnyPacket::size(); ++i, ++index)
                buffer[i] = T(calc(index));
            AnyPacket packet{};
            packet.load(buffer);
            return packet;
        }
    }

    template<class Derived>
    template<class AnyPacket>
    inline AnyPacket RValueVector<Derived>::packetPartial(size_t index, size_t count) const {
        using T = Traits<AnyPacket>::ScalarType;
        assert(index + count <= getLength() && "[Error]: Index out of range");
        if constexpr (Diffable<T>) {
            using ValuePacket = AnyPacket::ValueType;
            if constexpr (isForwardDiff) {
                using GradPacket = AnyPacket::GradType;
                const auto& x = Base::getDerived();
                auto values = toValueVector(x).template packetPartial<ValuePacket>(index, count);
                auto grads = toGradVector(x).template packetPartial<GradPacket>(index, count);
                return AnyPacket(std::move(values), std::move(grads));
            }
            else {
                using ValueType = T::ValueType;
                ValueType values[AnyPacket::size()];
                for (size_t i = 0; i < AnyPacket::size(); ++i, ++index)
                    values[i] = i < count ? ValueType(calc(index)) : ValueType(0);
                ValuePacket packet{};
                packet.load(values);
                return AnyPacket(std::move(packet));
            }
        }
        else {
            T buffer[AnyPacket::size()];
            for (size_t i = 0; i < AnyPacket::size(); ++i, ++index)
                buffer[i] = i < count ? T(calc(index)) : T(0);
            AnyPacket packet{};
            packet.load_partial(buffer, count);
            return packet;
        }
    }

    template<class Derived>
    inline auto RValueVector<Derived>::format() const {
        return FormatedVector<Derived>(Base::getDerived());
    }

    template<class Derived>
    inline auto RValueVector<Derived>::transpose() const noexcept {
        return TransposeVector<Derived>(Base::getDerived());
    }

    template<class Derived>
    inline auto RValueVector<Derived>::conjugate() const noexcept {
        return ConjugateVector<Derived>(Base::getDerived());
    }

    template<class Derived>
    inline auto RValueVector<Derived>::hermite() const noexcept {
        return HermiteVector<Derived>(Base::getDerived());
    }

    template<class Derived>
    inline RValueVector<Derived>::RealType RValueVector<Derived>::norm1() const {
        return abs(Base::getDerived()).sum();
    }

    template<class Derived>
    inline RValueVector<Derived>::RealType RValueVector<Derived>::norm2() const {
        return norm();
    }

    template<class Derived>
    inline RValueVector<Derived>::RealType RValueVector<Derived>::norm() const {
        return sqrt(Base::getDerived().squaredNorm());
    }

    template<class Derived>
    inline RValueVector<Derived>::RealType RValueVector<Derived>::squaredNorm() const {
        return toSquaredNormVector(Base::getDerived()).sum();
    }

    template<class Derived>
    RValueVector<Derived>::RealType RValueVector<Derived>::lnSquaredNorm() const {
        const auto& derived = Base::getDerived();
        const RealType maxabs = abs(derived).max();
        assert(maxabs > RealType(std::numeric_limits<ScalarType>::min()) && "[Error]: ln(0) is not allowed");
        const RealType factor = reciprocal(maxabs);
        return ln((derived * factor).squaredNorm()) + RealType(2) * ln(maxabs);
    }

    template<class Derived>
    inline RValueVector<Derived>::RealType RValueVector<Derived>::normInf() const {
        return abs(Base::getDerived()).max();
    }

    template<class Derived>
    RValueVector<Derived>::ScalarType RValueVector<Derived>::max() const {
        static_assert(!ScalarType::isComplex, "[Error]: Compare between complex number is ill defined");

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
    RValueVector<Derived>::ScalarType RValueVector<Derived>::min() const {
        static_assert(!ScalarType::isComplex, "[Error]: Compare between complex number is ill defined");

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
    auto RValueVector<Derived>::sum() const -> CoDiff<ScalarType> {
        assert(getLength() != 0 && "[Error]: Sum of a empty vector is not well defined");
        if constexpr (Internal::EnableSIMD<Derived>::value && !isReverseDiff) {
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
            co_return buffer.sum();
        }
        else if constexpr (isReverseDiff) {
            size_t i = 0;
            ValueType v = 0;
            auto elems = co_for([this, i]() { return i < getLength(); }, [&]() { ++i; }, [&, this, i]() {
                auto elem = calc(i);
                v += elem.value();
                return elem;
            });
            auto result = co_yield std::move(v);
            for (auto& elem : elems)
                elem.reverse(result.grad());
        }
        else {
            auto result = ScalarType(0);
            for(size_t i = 0; i < getLength(); ++i)
                result += calc(i);
            co_return std::move(result);
        }
    }

    template<class Derived>
    RValueVector<Derived>::ScalarType RValueVector<Derived>::lnSumExp() const {
        const Derived& v = Base::getDerived();
        ValueType m;
        if constexpr (isComplex)
            m = abs(toValueVector(v)).max();
        else
            m = max().value();
        return ln(exp(v - m).sum()) + m;
    }

    template<class Derived>
    RValueVector<Derived>::ScalarType RValueVector<Derived>::prod() const {
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
    template<Vector V>
    inline auto RValueVector<Derived>::crossProduct(const V& v) const noexcept {
        return CrossProduct<Derived, V>(Base::getDerived(), v);
    }

    template<class Derived>
    template<Vector V>
    RValueVector<Derived>::ScalarType RValueVector<Derived>::angleTo(const V& v) const noexcept {
        return arccos(Base::getDerived() * v / (norm() * v.norm()));
    }

    template<class Derived>
    template<size_t Length>
    inline auto RValueVector<Derived>::head(size_t to) noexcept {
        return BlockType<Length>(Base::getDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    inline const auto RValueVector<Derived>::head(size_t to) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    inline auto RValueVector<Derived>::tail(size_t from) noexcept {
        return BlockType<Length>(Base::getDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    inline const auto RValueVector<Derived>::tail(size_t from) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    inline auto RValueVector<Derived>::segment(size_t from, size_t to) noexcept {
        return BlockType<Length>(Base::getDerived(), from, to);
    }

    template<class Derived>
    template<size_t Length>
    inline const auto RValueVector<Derived>::segment(size_t from, size_t to) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), from, to);
    }

    template<class Derived>
    inline auto RValueVector<Derived>::reverse() noexcept {
        return ReverseVector<Derived>(Base::getDerived());
    }

    template<class Derived>
    inline const auto RValueVector<Derived>::reverse() const noexcept {
        return ReverseVector<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueVector<Derived>::reals() const noexcept {
        return RealVector<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueVector<Derived>::imags() const noexcept {
        return ImagVector<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueVector<Derived>::squaredNorms() const noexcept {
        return SquaredNormVector<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueVector<Derived>::norms() const noexcept {
        return NormVector<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueVector<Derived>::values() const noexcept {
        return ValueVector<Derived>(Base::getDerived());
    }

    template<class Derived>
    template<int GradOrder>
    auto RValueVector<Derived>::grads() const noexcept {
        return GradVector<Derived, GradOrder>(Base::getDerived());
    }

    template<Vector T1, Vector T2>
    bool vectorNear(const T1& v1, const T2& v2, double precision) {
        using ScalarType = Internal::BinaryScalarOpRtnTy<typename T1::ScalarType, typename T2::ScalarType>::Type;
        assert(v1.getLength() == v2.getLength());
        for (size_t i = 0; i < v1.getLength(); ++i)
            if (!scalarNear(ScalarType(v1.calc(i)), ScalarType(v2.calc(i)), precision))
                return false;
        return true;
    }
}
