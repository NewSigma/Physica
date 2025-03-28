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

namespace Physica {
    template<class Derived>
    template<Vector V, class Executor>
    inline void RValueVector<Derived>::assign(V& v) const {
        constexpr static size_t Size1 = SizeAtCompile;
        constexpr static size_t Size2 = V::SizeAtCompile;
        assign_check(v);

        assert(getLength() == v.getLength() && "[Error]: Size mismatch between two vector");
        if constexpr (Internal::EnableSIMD<Derived, V>::value && !isReverseDiff) {
            constexpr static size_t Size = Size1 > Size2 ? Size1 : Size2;
            assign_simd<V, Executor, Size>(v);
        }
        else
            assign_for<V, Executor>(v);
    }

    template<class Derived>
    template<Vector V, class Executor>
    inline void RValueVector<Derived>::assign_add(V& v) const {
        constexpr static size_t Size1 = SizeAtCompile;
        constexpr static size_t Size2 = V::SizeAtCompile;
        assign_check(v);

        assert(getLength() == v.getLength() && "[Error]: Size mismatch between two vector");
        if constexpr (Internal::EnableSIMD<Derived, V>::value && !isReverseDiff) {
            constexpr static size_t Size = Size1 > Size2 ? Size1 : Size2;
            assign_add_simd<V, Size>(v);
        }
        else
            assign_add_for<V, Executor>(v);
    }

    template<class Derived>
    template<Packet Pack>
    inline Pack RValueVector<Derived>::packet(size_t index) const {
        using U = Traits<Pack>::ScalarType;
        assert(index + Pack::size() <= getLength() && "[Error]: Index out of range");
        if constexpr (Diffable<U>) {
            using ValuePacket = Pack::ValueType;
            if constexpr (isForwardDiff) {
                using GradPacket = Pack::GradType;
                const auto& x = Base::getDerived();
                auto values = x.values().template packet<ValuePacket>(index);
                auto grads = x.grads().template packet<GradPacket>(index);
                return Pack(std::move(values), std::move(grads));
            }
            else {
                using Uv = U::ValueType;
                Uv values[Pack::size()];
                for (size_t i = 0; i < Pack::size(); ++i, ++index)
                    values[i] = Uv(calc(index));
                ValuePacket packet{};
                packet.load(values);
                return Pack(std::move(packet));
            }
        }
        else {
            U buffer[Pack::size()];
            for (size_t i = 0; i < Pack::size(); ++i, ++index)
                buffer[i] = U(calc(index));
            Pack packet{};
            packet.load(buffer);
            return packet;
        }
    }

    template<class Derived>
    template<Packet Pack>
    inline Pack RValueVector<Derived>::packetPartial(size_t index, size_t count) const {
        using U = Traits<Pack>::ScalarType;
        assert(index + count <= getLength() && "[Error]: Index out of range");
        if constexpr (Diffable<U>) {
            using ValuePacket = Pack::ValueType;
            if constexpr (isForwardDiff) {
                using GradPacket = Pack::GradType;
                const auto& x = Base::getDerived();
                auto values = x.values().template packetPartial<ValuePacket>(index, count);
                auto grads = x.grads().template packetPartial<GradPacket>(index, count);
                return Pack(std::move(values), std::move(grads));
            }
            else {
                using Uv = U::ValueType;
                Uv values[Pack::size()];
                for (size_t i = 0; i < Pack::size(); ++i, ++index)
                    values[i] = i < count ? Uv(calc(index)) : Uv(0);
                ValuePacket packet{};
                packet.load(values);
                return Pack(std::move(packet));
            }
        }
        else {
            U buffer[Pack::size()];
            for (size_t i = 0; i < Pack::size(); ++i, ++index)
                buffer[i] = i < count ? U(calc(index)) : U(0);
            Pack packet{};
            packet.load_partial(buffer, count);
            return packet;
        }
    }

    template<class Derived>
    template<Vector V1, Vector V2>
    void RValueVector<Derived>::reverse(const V1&, const V2& grad) const noexcept requires(isReverseDiff) {
        Base::getDerived().reverse(grad);
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
    inline auto RValueVector<Derived>::norm1() const -> CoDiff<Tr> {
        return abs(Base::getDerived()).sum();
    }

    template<class Derived>
    inline auto RValueVector<Derived>::norm2() const -> CoDiff<Tr> {
        return norm();
    }

    template<class Derived>
    inline auto RValueVector<Derived>::norm() const -> CoDiff<Tr> {
        return sqrt(squaredNorm());
    }

    template<class Derived>
    inline auto RValueVector<Derived>::squaredNorm() const -> CoDiff<Tr> {
        return squaredNorms().sum();
    }

    template<class Derived>
    auto RValueVector<Derived>::lnSquaredNorm() const -> Tr {
        const auto& derived = Base::getDerived();
        const Tr maxabs = abs(derived).max();
        assert(maxabs > Trv(std::numeric_limits<T>::min()) && "[Error]: ln(0) is not allowed");
        const Tr factor = reciprocal(maxabs);
        return ln((derived * factor).squaredNorm()) + Tr(2) * ln(maxabs);
    }

    template<class Derived>
    inline auto RValueVector<Derived>::normInf() const -> Tr {
        return abs(Base::getDerived()).max();
    }

    template<class Derived>
    auto RValueVector<Derived>::max() const -> CoDiff<T> {
        static_assert(!T::isComplex, "[Error]: Compare between complex number is ill defined");
        assert(getLength() != 0);

        constexpr bool EnableSIMD = Internal::EnableSIMD<Derived>::value;
        if constexpr (isReverseDiff) {
            size_t index = 0;
            Tv elem = calc_value(0), temp;
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
            T result = calc(0);
            for(size_t i = 1; i < getLength(); ++i) {
                T temp = calc(i);
                if (result < temp)
                    result = temp;
            }
            co_return std::move(result);
        }
        else {
            const auto& v = Base::getDerived();
            PacketType buffer(std::numeric_limits<T>::lowest());
            T result;
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
                        result = std::max<T>(result, v.calc(i + j));
                }
            }
            co_return std::move(result);
        }
    }

    template<class Derived>
    auto RValueVector<Derived>::min() const -> CoDiff<T> {
        static_assert(!T::isComplex, "[Error]: Compare between complex number is ill defined");
        assert(getLength() != 0);

        constexpr bool EnableSIMD = Internal::EnableSIMD<Derived>::value;
        if constexpr (isReverseDiff) {
            size_t index = 0;
            Tv elem = calc_value(0), temp;
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
            PacketType buffer(std::numeric_limits<T>::max());
            T result;
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
                        result = std::min<T>(result, v.calc(i + j));
                }
            }
            co_return std::move(result);
        }
        else {
            T result = calc(0);
            for(size_t i = 1; i < getLength(); ++i) {
                T temp = calc(i);
                if (result > temp)
                    result = temp;
            }
            co_return std::move(result);
        }
    }

    template<class Derived>
    auto RValueVector<Derived>::sum() const -> CoDiff<T> {
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
            auto result = T(0);
            for(size_t i = 0; i < getLength(); ++i)
                result += calc(i);
            co_return std::move(result);
        }
    }

    template<class Derived>
    auto RValueVector<Derived>::mean() const -> CoDiff<T> {
        return sum() / Trv(getLength());
    }

    template<class Derived>
    auto RValueVector<Derived>::variance() const -> CoDiff<T> {
        const auto& x = Base::getDerived();
        const size_t length = getLength();
        const auto x_mean = mean();
        const auto expr = x - x_mean;
        const auto expr2 = square(expr);
        auto result = expr2.sum() / Trv(length);
        if constexpr (isReverseDiff) {
            auto tmp = co_yield result.value();
            result.reverse(tmp.grad());
        }
        else
            co_return std::move(result);
    }

    template<class Derived>
    auto RValueVector<Derived>::variance(const T& prior_mean) const -> T {
        const auto& x = Base::getDerived();
        const size_t length = getLength();
        return (x - prior_mean).squaredNorm() / Trv(length - 1);
    }

    template<class Derived>
    auto RValueVector<Derived>::deviation() const -> T {
        return sqrt(variance());
    }

    template<class Derived>
    auto RValueVector<Derived>::deviation(const T& prior_mean) const -> T {
        return sqrt(variance(prior_mean));
    }

    template<class Derived>
    auto RValueVector<Derived>::lnSumExp() const -> CoDiff<T> {
        const Derived& v = Base::getDerived();
        Tv m;
        if constexpr (isComplex)
            m = values().reals().max();
        else
            m = values().max();

        auto expr1 = v - m;
        auto expr2 = exp(expr1);
        auto y = ln(expr2.sum() + Trv(std::numeric_limits<Trv>::min())) + m; // Add min() to avoid ln(0)
        if constexpr (isReverseDiff) {
            auto tmp = co_yield y.value();
            y.reverse_final(tmp.grad());
        }
        else
            co_return std::move(y);
    }

    template<class Derived>
    auto RValueVector<Derived>::crossEntropy(size_t index) const -> CoDiff<T> {
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
    auto RValueVector<Derived>::lnSoftmax(size_t index) const -> CoDiff<T> {
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
    auto RValueVector<Derived>::softmax(size_t index) const -> CoDiff<T> {
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
    auto RValueVector<Derived>::prod() const -> CoDiff<T> {
        assert(getLength() != 0);
        if constexpr (isReverseDiff) {
            Tv result = calc_value(0);
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
    auto RValueVector<Derived>::angleTo(const V& v) const noexcept -> T {
        return arccos(Base::getDerived() * v / (norm() * v.norm()));
    }

    template<class Derived>
    template<size_t Length>
    inline auto RValueVector<Derived>::head(size_t to) & noexcept {
        return BlockType<Length>(Base::getDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    inline const auto RValueVector<Derived>::head(size_t to) const& noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    inline auto RValueVector<Derived>::tail(size_t from) & noexcept {
        return BlockType<Length>(Base::getDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    inline const auto RValueVector<Derived>::tail(size_t from) const& noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    inline auto RValueVector<Derived>::segment(size_t from, size_t to) & noexcept {
        return BlockType<Length>(Base::getDerived(), from, to);
    }

    template<class Derived>
    template<size_t Length>
    inline const auto RValueVector<Derived>::segment(size_t from, size_t to) const& noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), from, to);
    }

    template<class Derived>
    auto RValueVector<Derived>::reals() const noexcept -> RealsRtnTy {
        return RealsRtnTy(Base::getDerived());
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
    auto RValueVector<Derived>::values() const noexcept -> ValuesRtnTy {
        return ValuesRtnTy(Base::getDerived());
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
    template<Vector V>
    __host__ __device__ void RValueVector<Derived>::assign_check(const V&) noexcept {
        constexpr static size_t Size1 = SizeAtCompile;
        constexpr static size_t Size2 = V::SizeAtCompile;
        static_assert(Size1 == Dynamic || Size2 == Dynamic || Size1 == Size2, "[Error]: Size mismatch between two vector");
        static_assert(V::isComplex || !isComplex, "[Error]: Cannot convert a complex to a real");
        static_assert(Diffable<V> || !Diffable<This>, "[Error]: Assign a diffable vector to normal vector discards grads");
    }

    template<class Derived>
    template<int GradOrder>
    auto RValueVector<Derived>::grads_impl() const noexcept {
        return GradVector<Derived, GradOrder>(Base::getDerived());
    }

    template<class Derived>
    template<Vector V, class Executor>
    inline void RValueVector<Derived>::assign_for(V& v) const {
        Executor::parallel_for([&, this](size_t i) {
            if constexpr (isReverseDiff)
                v[i] = calc_value(i);
            else
                v[i] = calc(i);
        }, getLength(), Executor::getNumThread()).wait();
    }

    template<class Derived>
    template<Vector V, class Executor, size_t Size>
    inline void RValueVector<Derived>::assign_simd(V& v) const {
        using Pack = BestPacket<typename V::ScalarType, Size>::Type;
        constexpr static size_t PacketSize = Pack::size();
        const auto& v0 = Base::getDerived();
        if constexpr (Size != Dynamic) {
            constexpr size_t to = Size / PacketSize * PacketSize;
            for (size_t i = 0; i < to; i += PacketSize)
                v.writePacket(i, v0.template packet<Pack>(i));
            
            constexpr size_t i = Size - Size % PacketSize;
            if constexpr (i != Size) {
                constexpr size_t count = Size - i;
                v.writePacketPartial(i, count, v0.template packetPartial<Pack>(i, count));
            }
        }
        else {
            const size_t length = getLength();
            if (length == 0)
                return;

            const size_t to = length / PacketSize * PacketSize;
            if constexpr (std::is_same<Executor, SeqExecutor>::value) {
                size_t i = 0;
                for (; i < to; i += PacketSize)
                    v.writePacket(i, v0.template packet<Pack>(i));

                if (to != length) {
                    const size_t count = length - i;
                    v.writePacketPartial(i, count, v0.template packetPartial<Pack>(i, count));
                }
            }
            else {
                const size_t numLoop = to / PacketSize;
                auto future = Executor::parallel_for([&, this](size_t i) {
                    const size_t i1 = i * PacketSize;
                    v.writePacket(i1, v0.template packet<Pack>(i1));
                }, numLoop, Executor::getNumThread());

                if (to != length) {
                    const size_t count = length % PacketSize;
                    v.writePacketPartial(to, count, v0.template packetPartial<Pack>(to, count));
                }
                future.wait();
            }
        }
    }

    template<class Derived>
    template<Vector V, class Executor>
    inline void RValueVector<Derived>::assign_add_for(V& v) const {
        Executor::parallel_for([&, this](size_t i) {
            if constexpr (isReverseDiff)
                v[i] += calc_value(i);
            else
                v[i] += calc(i);
        }, getLength(), Executor::getNumThread()).wait();
    }

    template<class Derived>
    template<Vector V, size_t Size>
    inline void RValueVector<Derived>::assign_add_simd(V& v) const {
        using Pack = BestPacket<typename V::ScalarType, Size>::Type;
        constexpr static size_t PacketSize = Pack::size();
        const auto& v0 = Base::getDerived();
        if constexpr (Size != Dynamic) {
            constexpr size_t to = Size / PacketSize * PacketSize;
            for (size_t i = 0; i < to; i += PacketSize) {
                const Pack sum = v.template packet<Pack>(i) + v0.template packet<Pack>(i);
                v.writePacket(i, sum);
            }

            constexpr size_t i = Size - Size % PacketSize;
            if constexpr (i != Size) {
                constexpr size_t count = Size - i;
                const Pack sum = v.template packetPartial<Pack>(i, count) + v0.template packetPartial<Pack>(i, count);
                v.writePacketPartial(i, count, sum);
            }
        }
        else {
            const size_t length = v.getLength();
            size_t i = 0;
            const size_t to = length / PacketSize * PacketSize;
            for (; i < to; i += PacketSize) {
                const Pack sum = v.template packet<Pack>(i) + v0.template packet<Pack>(i);
                v.writePacket(i, sum);
            }
            if (to != length) {
                const size_t count = length - i;
                const Pack sum = v.template packetPartial<Pack>(i, count) + v0.template packetPartial<Pack>(i, count);
                v.writePacketPartial(i, count, sum);
            }
        }
    }

    template<Vector T1, Vector T2>
    bool vectorNear(const T1& v1, const T2& v2, double precision) {
        using T = Internal::BinaryScalarOpRtnTy<typename T1::ScalarType, typename T2::ScalarType>::Type;
        assert(v1.getLength() == v2.getLength());
        for (size_t i = 0; i < v1.getLength(); ++i)
            if (!scalarNear(T(v1.calc(i)), T(v2.calc(i)), precision))
                return false;
        return true;
    }
}
