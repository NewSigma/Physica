/*
 * Copyright 2022-2025 Weibo He.
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
        template<Vector T1, Vector T2, class Executor>
        class AssignImpl {
            constexpr static size_t Size1 = T1::SizeAtCompile;
            constexpr static size_t Size2 = T2::SizeAtCompile;
            constexpr static size_t SizeAtCompile = Size1 > Size2 ? Size1 : Size2;
            static_assert(Size1 == Dynamic || Size2 == Dynamic || Size1 == Size2, "[Error]: Size mismatch between two vector");
            static_assert(T1::isComplex || !T2::isComplex, "[Error]: Cannot convert a complex to a real");
            static_assert(Diffable<T1> || !Diffable<T2>, "[Error]: Assign a diffable vector to normal vector discards grads");
        public:
            using ScalarType = T1::ScalarType;
            using AnyPacket = BestPacket<ScalarType, SizeAtCompile>::Type;
            constexpr static size_t PacketSize = AnyPacket::size();
            constexpr static bool isReverseDiff = ReverseDiff<ScalarType>;

            inline static void run(T1& v1, const T2& v2) {
                assert(v1.getLength() == v2.getLength() && "[Error]: Size mismatch between two vector");
                if constexpr (Internal::EnableSIMD<T1, T2>::value && !isReverseDiff)
                    run_simd(v1, v2);
                else
                    run_for(v1, v2);
            }
        private:
            inline static void run_for(T1& v1, const T2& v2) {
                Executor::parallel_for([&](size_t i) {
                    if constexpr (isReverseDiff)
                        v1[i] = v2.calc_value(i);
                    else
                        v1[i] = v2.calc(i);
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
                    if constexpr (std::is_same<Executor, SeqExecutor>::value) {
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
    template<Vector V, class Executor>
    inline void RValueVector<Derived>::assign(V& v) const {
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
                auto values = x.values().template packet<ValuePacket>(index);
                auto grads = x.grads().template packet<GradPacket>(index);
                return AnyPacket(std::move(values), std::move(grads));
            }
            else {
                using Tv = T::ValueType;
                Tv values[AnyPacket::size()];
                for (size_t i = 0; i < AnyPacket::size(); ++i, ++index)
                    values[i] = Tv(calc(index));
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
                auto values = x.values().template packetPartial<ValuePacket>(index, count);
                auto grads = x.grads().template packetPartial<GradPacket>(index, count);
                return AnyPacket(std::move(values), std::move(grads));
            }
            else {
                using Tv = T::ValueType;
                Tv values[AnyPacket::size()];
                for (size_t i = 0; i < AnyPacket::size(); ++i, ++index)
                    values[i] = i < count ? Tv(calc(index)) : Tv(0);
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
    inline auto RValueVector<Derived>::norm1() const -> CoDiff<RealType> {
        return abs(Base::getDerived()).sum();
    }

    template<class Derived>
    inline auto RValueVector<Derived>::norm2() const -> CoDiff<RealType> {
        return norm();
    }

    template<class Derived>
    inline auto RValueVector<Derived>::norm() const -> CoDiff<RealType> {
        return sqrt(squaredNorm());
    }

    template<class Derived>
    inline auto RValueVector<Derived>::squaredNorm() const -> CoDiff<RealType> {
        return squaredNorms().sum();
    }

    template<class Derived>
    auto RValueVector<Derived>::lnSquaredNorm() const -> RealType {
        const auto& derived = Base::getDerived();
        const RealType maxabs = abs(derived).max();
        assert(maxabs > RealType(std::numeric_limits<ScalarType>::min()) && "[Error]: ln(0) is not allowed");
        const RealType factor = reciprocal(maxabs);
        return ln((derived * factor).squaredNorm()) + RealType(2) * ln(maxabs);
    }

    template<class Derived>
    inline auto RValueVector<Derived>::normInf() const -> RealType {
        return abs(Base::getDerived()).max();
    }

    template<class Derived>
    auto RValueVector<Derived>::max() const -> CoDiff<ScalarType> {
        static_assert(!ScalarType::isComplex, "[Error]: Compare between complex number is ill defined");
        assert(getLength() != 0);

        constexpr bool EnableSIMD = Internal::EnableSIMD<Derived>::value;
        if constexpr (isReverseDiff) {
            size_t index = 0;
            ValueType elem = calc_value(0), temp;
            for(size_t i = 1; i < getLength(); ++i) {
                temp = calc_value(i);
                if (elem < temp) {
                    elem = temp;
                    index = i;
                }
            }
            const auto result = co_yield calc_value(index);
            calc(index).reverse(result.grad());
        }
        else if constexpr (!EnableSIMD) {
            ScalarType result = calc(0);
            for(size_t i = 1; i < getLength(); ++i) {
                ScalarType temp = calc(i);
                if (result < temp)
                    result = temp;
            }
            co_return std::move(result);
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
                        result = std::max<ScalarType>(result, v.calc(i + j));
                }
            }
            co_return std::move(result);
        }
    }

    template<class Derived>
    auto RValueVector<Derived>::min() const -> CoDiff<ScalarType> {
        static_assert(!ScalarType::isComplex, "[Error]: Compare between complex number is ill defined");
        assert(getLength() != 0);

        constexpr bool EnableSIMD = Internal::EnableSIMD<Derived>::value;
        if constexpr (isReverseDiff) {
            size_t index = 0;
            ValueType elem = calc_value(0), temp;
            for(size_t i = 1; i < getLength(); ++i) {
                temp = calc_value(i);
                if (elem > temp) {
                    elem = temp;
                    index = i;
                }
            }
            const auto result = co_yield calc_value(index);
            calc(index).reverse(result.grad());
        }
        else if constexpr (EnableSIMD) {
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
                        result = std::min<ScalarType>(result, v.calc(i + j));
                }
            }
            co_return std::move(result);
        }
        else {
            ScalarType result = calc(0);
            for(size_t i = 1; i < getLength(); ++i) {
                ScalarType temp = calc(i);
                if (result > temp)
                    result = temp;
            }
            co_return std::move(result);
        }
    }

    template<class Derived>
    auto RValueVector<Derived>::sum() const -> CoDiff<ScalarType> {
        assert(getLength() != 0 && "[Error]: Sum of a empty vector is not well defined");
        if constexpr (isReverseDiff) {
            auto result = co_yield values().sum();
            Base::getDerived().reverse(result.grad());
        }
        else if constexpr (Internal::EnableSIMD<Derived>::value) {
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
        else {
            auto result = ScalarType(0);
            for(size_t i = 0; i < getLength(); ++i)
                result += calc(i);
            co_return std::move(result);
        }
    }

    template<class Derived>
    auto RValueVector<Derived>::lnSumExp() const -> CoDiff<ScalarType> {
        const Derived& v = Base::getDerived();
        ValueType m;
        if constexpr (isComplex)
            m = values().reals().max();
        else
            m = values().max();

        auto expr1 = v - m;
        auto expr2 = exp(expr1);
        auto y = ln(expr2.sum() + ValueType(std::numeric_limits<ValueType>::min())) + m; // Add min() to avoid ln(0)
        if constexpr (isReverseDiff) {
            auto tmp = co_yield y.value();
            y.reverse_final(tmp.grad());
        }
        else
            co_return std::move(y);
    }

    template<class Derived>
    auto RValueVector<Derived>::crossEntropy(size_t index) const -> CoDiff<ScalarType> {
        assert(index < getLength() && "[Error]: Index overflow");
        const auto& v = Base::getDerived();
        const auto vi = calc(index);
        const auto delta = v - vi;
        auto y = delta.lnSumExp();
        if constexpr (isReverseDiff) {
            auto tmp = co_yield y.value();
            y.reverse_final(tmp.grad());
        }
        else
            co_return std::move(y);
    }

    template<class Derived>
    auto RValueVector<Derived>::lnSoftmax(size_t index) const -> CoDiff<ScalarType> {
        assert(index < getLength() && "[Error]: Index overflow");
        auto y = -crossEntropy(index);
        if constexpr (isReverseDiff) {
            auto tmp = co_yield y.value();
            y.reverse(tmp.grad());
        }
        else
            co_return std::move(y);
    }

    template<class Derived>
    auto RValueVector<Derived>::softmax(size_t index) const -> CoDiff<ScalarType> {
        assert(index < getLength() && "[Error]: Index overflow");
        auto y = exp(lnSoftmax(index));
        if constexpr (isReverseDiff) {
            auto tmp = co_yield y.value();
            y.reverse_final(tmp.grad());
        }
        else
            co_return std::move(y);
    }

    template<class Derived>
    auto RValueVector<Derived>::prod() const -> CoDiff<ScalarType> {
        assert(getLength() != 0);
        if constexpr (isReverseDiff) {
            ValueType result = calc_value(0);
            for(size_t i = 1; i < getLength(); ++i)
                result *= calc_value(i);

            auto tmp = co_yield std::move(result);
            const auto& v = Base::getDerived();
            v.reverse(reciprocal(v.values()) * (tmp.value() * tmp.grad()));
        }
        else {
            auto result = calc(0);
            for(size_t i = 1; i < getLength(); ++i)
                result *= calc(i);
            co_return std::move(result);
        }
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
    auto RValueVector<Derived>::angleTo(const V& v) const noexcept -> ScalarType {
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
        if constexpr (isReverseDiff)
            return Base::getDerived().template grads<GradOrder>();
        else
            return grads_impl<GradOrder>();
    }

    template<class Derived>
    template<int GradOrder>
    auto RValueVector<Derived>::grads_impl() const noexcept {
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
