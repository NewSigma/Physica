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
    template<ExecutePolicy P>
    void RValueVector<Derived>::assign(Vector auto&& v) const noexcept {
        using V = std::remove_cvref<decltype(v)>::type;
        constexpr static size_t Length1 = SizeAtCompile;
        constexpr static size_t Length2 = V::SizeAtCompile;
        assert_assign(v);

        assert(getLength() == v.getLength() && "[Error]: Length mismatch between two vector");
        if constexpr (Internal::EnableSIMD<Derived, V>::value && !isReverseDiff) {
            constexpr static size_t Length = Length1 > Length2 ? Length1 : Length2;
            assign_simd<V, P, Length>(v);
        }
        else
            assign_for<V, P>(v);
    }

    template<class Derived>
    template<ExecutePolicy P>
    void RValueVector<Derived>::assign_base(Vector auto&& v) const noexcept {
        assign<P>(v);
    }

    template<class Derived>
    template<ExecutePolicy P>
    void RValueVector<Derived>::assign_add(Vector auto& v) const noexcept {
        using V = std::remove_cvref<decltype(v)>::type;
        constexpr static size_t Length1 = SizeAtCompile;
        constexpr static size_t Length2 = V::SizeAtCompile;
        assert_assign(v);

        assert(getLength() == v.getLength() && "[Error]: Length mismatch between two vector");
        if constexpr (Internal::EnableSIMD<Derived, V>::value && !isReverseDiff) {
            constexpr static size_t Length = Length1 > Length2 ? Length1 : Length2;
            assign_add_simd<V, Length>(v);
        }
        else
            assign_add_for<V, P>(v);
    }

    template<class Derived>
    void RValueVector<Derived>::assert_assign(const Vector auto& target) const noexcept {
        static_assert_assign(target);

        using V = std::remove_cvref<decltype(target)>::type;
        constexpr size_t Size1 = SizeAtCompile;
        constexpr size_t Size2 = V::SizeAtCompile;
        if constexpr (Size1 == Dynamic && Size2 == Dynamic)
            assert(getLength() == target.getLength() && "[Error]: Size mismatch between two vector");
    }

    template<class Derived>
    void RValueVector<Derived>::assert_assign_mkl(const Vector auto& target) const noexcept {
        assert_assign(target);
        static_assert(Internal::EnableMKL<Derived, decltype(target)>::value, "[Error]: Cannot apply MKL to this expr");
    }

    template<class Derived>
    decltype(auto) RValueVector<Derived>::calc(size_t index) const noexcept {
        return Base::getDerived().calc(index);
    }

    template<class Derived>
    decltype(auto) RValueVector<Derived>::calc_value(size_t index) const noexcept {
        return Base::getDerived().calc_value(index);
    }

    template<class Derived>
    template<Packet Pack>
    Pack RValueVector<Derived>::packet(size_t index) const noexcept {
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
                Array<Uv, Pack::size()> values{};
                for (size_t i = 0; i < Pack::size(); ++i, ++index)
                    values[i] = Uv(calc(index));
                ValuePacket packet{};
                packet.load(values.data());
                return Pack(std::move(packet));
            }
        }
        else {
            Array<U, Pack::size()> buffer{};
            for (size_t i = 0; i < Pack::size(); ++i, ++index)
                buffer[i] = U(calc(index));
            Pack packet{};
            packet.load(buffer.data());
            return packet;
        }
    }

    template<class Derived>
    template<Packet Pack>
    Pack RValueVector<Derived>::packetPartial(size_t index, size_t count) const noexcept {
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
                Array<Uv, Pack::size()> values{};
                for (size_t i = 0; i < Pack::size(); ++i, ++index)
                    values[i] = i < count ? Uv(calc(index)) : Uv(0);
                ValuePacket packet{};
                packet.load(values.data());
                return Pack(std::move(packet));
            }
        }
        else {
            Array<U, Pack::size()> buffer{};
            for (size_t i = 0; i < Pack::size(); ++i, ++index)
                buffer[i] = i < count ? U(calc(index)) : U(0);
            Pack packet{};
            packet.load_partial(buffer.data(), count);
            return packet;
        }
    }

    template<class Derived>
    void RValueVector<Derived>::reverse(const Vector auto&, const Vector auto& grad) const noexcept requires(isReverseDiff) {
        Base::getDerived().reverse(grad);
    }

    template<class Derived>
    template<size_t Length>
    auto RValueVector<Derived>::head(size_t to) & noexcept {
        return BlockType<Length>(Base::getDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    const auto RValueVector<Derived>::head(size_t to) const& noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    auto RValueVector<Derived>::tail(size_t from) & noexcept {
        return BlockType<Length>(Base::getDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    const auto RValueVector<Derived>::tail(size_t from) const& noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    auto RValueVector<Derived>::segment(size_t from, size_t to) & noexcept {
        return BlockType<Length>(Base::getDerived(), from, to);
    }

    template<class Derived>
    template<size_t Length>
    const auto RValueVector<Derived>::segment(size_t from, size_t to) const& noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), from, to);
    }

    template<class Derived>
    auto RValueVector<Derived>::format() const {
        return FormatedVector<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueVector<Derived>::transpose() const noexcept {
        return TransposeVector<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueVector<Derived>::conjugate() const noexcept -> ConjugateRtnTy {
        if constexpr (isComplex)
            return Conjugate<Derived>(Base::getDerived());
        else
            return Base::getDerived();
    }

    template<class Derived>
    auto RValueVector<Derived>::hermite() const noexcept {
        return HermiteVector<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueVector<Derived>::norm1() const noexcept -> CoDiff<Tr> {
        return abs(Base::getDerived()).sum();
    }

    template<class Derived>
    auto RValueVector<Derived>::norm2() const noexcept -> CoDiff<Tr> {
        return sqrt(squaredNorm());
    }

    template<class Derived>
    auto RValueVector<Derived>::norm() const noexcept {
        return Base::getDerived().norm2();
    }

    template<class Derived>
    auto RValueVector<Derived>::squaredNorm() const noexcept -> CoDiff<Tr> {
        return squaredNorms().sum();
    }

    template<class Derived>
    auto RValueVector<Derived>::lnSquaredNorm() const -> Tr {
        const auto& derived = Base::getDerived();
        const Tr maxabs = abs(derived).max();
        // We require a small threhold to avoid latter ill-conditioned reciprocal.
        assert(maxabs > Trv(std::numeric_limits<T>::min()) && "[Error]: Vectors near zero are invalid");
        const Tr factor = reciprocal(maxabs);
        return ln((derived * factor).squaredNorm()) + Tr(2) * ln(maxabs);
    }

    template<class Derived>
    auto RValueVector<Derived>::normInf() const -> Tr {
        return abs(Base::getDerived()).max();
    }

    template<class Derived>
    auto RValueVector<Derived>::max() const noexcept -> CoDiff<T> {
        static_assert(!T::isComplex, "[Error]: Compare between complex number is ill defined");
        assert(getLength() != 0);

        constexpr bool EnableSIMD = Internal::EnableSIMD<Derived>::value;
        if constexpr (isReverseDiff) {
            size_t index = 0;
            Tv elem = calc_value(0);
            for(size_t i = 1; i < getLength(); ++i) {
                const Tv temp = calc_value(i);
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
    auto RValueVector<Derived>::min() const noexcept -> CoDiff<T> {
        static_assert(!T::isComplex, "[Error]: Compare between complex number is ill defined");
        assert(getLength() != 0);

        constexpr bool EnableSIMD = Internal::EnableSIMD<Derived>::value;
        if constexpr (isReverseDiff) {
            size_t index = 0;
            Tv elem = calc_value(0);
            for(size_t i = 1; i < getLength(); ++i) {
                const Tv temp = calc_value(i);
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
    auto RValueVector<Derived>::sum() const noexcept -> CoDiff<T> {
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
    auto RValueVector<Derived>::mean() const noexcept -> CoDiff<T> {
        return sum() / Trv(getLength());
    }

    template<class Derived>
    auto RValueVector<Derived>::mean_stable() const noexcept -> CoDiff<T> {
        auto result = T(0);
        const auto& v = Base::getDerived();
        for(size_t i = 0; i < getLength(); ++i)
           result.toNextMean(i, v.calc(i));

        if constexpr (isReverseDiff) {
            auto y = co_yield std::move(result);
            v.reverse(y.grad());
        }
        else
            co_return std::move(result);
    }

    template<class Derived>
    auto RValueVector<Derived>::variance() const noexcept -> CoDiff<T> {
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
    /**
     * More stable for large dataset. Prior version is not provided because they behave similarly at large dataset.
     */
    template<class Derived>
    auto RValueVector<Derived>::variance_stable() const -> T {
        using ScalarType = T::ScalarType;
        ScalarType result = 0;
        ScalarType mean = 0;
        for (size_t i = 0; i < getLength(); ++i)
            result.toNextVariance(mean, i, calc(i));
        return result;
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
    auto RValueVector<Derived>::deviation_stable() const -> T {
        return sqrt(variance_stable());
    }

    template<class Derived>
    auto RValueVector<Derived>::lnSumExp() const noexcept -> CoDiff<T> {
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
    auto RValueVector<Derived>::crossEntropy(size_t index) const noexcept -> CoDiff<T> {
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
    auto RValueVector<Derived>::lnSoftmax(size_t index) const noexcept -> CoDiff<T> {
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
    auto RValueVector<Derived>::softmax(size_t index) const noexcept -> CoDiff<T> {
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
    auto RValueVector<Derived>::prod() const noexcept -> CoDiff<T> {
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
    bool RValueVector<Derived>::isFinite() const {
        for(size_t i = 0; i < getLength(); ++i)
            if (!calc(i).isFinite())
                return false;
        return true;
    }

    template<class Derived>
    auto RValueVector<Derived>::crossProduct(const Vector auto& v) const noexcept {
        using V = std::remove_cvref_t<decltype(v)>;
        return CrossProduct<Derived, V>(Base::getDerived(), v);
    }

    template<class Derived>
    auto RValueVector<Derived>::angleTo(const Vector auto& v) const noexcept -> T {
        return arccos(Base::getDerived() * v / (norm() * v.norm()));
    }
    /**
     * The first element of \param target will be the factor to construct houseHolder matrix.
     * The other parts of \param target will be essential HouseHolder vector.
     * 
     * \return Norm of original vector
     * 
     * References:
     * [1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013:236-244
     * [2] Eigen; https://eigen.tuxfamily.org
     */
    template<class Derived>
    auto RValueVector<Derived>::householder(Vector auto& target) const -> Tr {
        using V = std::remove_cvref<decltype(target)>::type;
        constexpr size_t Length = std::max(SizeAtCompile, V::SizeAtCompile);
        constexpr size_t TailLength = Length > 0 ? (Length - 1) : Dynamic;
        assert(getLength() == target.getLength());
        assert(getLength() > 1 && "[Error]: Unnecessary householder call");

        const T v0 = calc(0);
        const Tr sourceNorm0 = v0.squaredNorm();
        const Tr squaredTailNorm = tail<TailLength>(1).squaredNorm();
        if (squaredTailNorm > std::numeric_limits<T>::min()) [[likely]] {
            const Tr norm = sqrt(squaredTailNorm + sourceNorm0);
            target[0] = Tr(1) + abs(v0) / norm;
            target.template tail<TailLength>(1) = tail<TailLength>(1) * reciprocal(v0 + unit(v0.value()) * norm);
            return norm;
        }
        else {
            const bool isZeroVector = sourceNorm0 < std::numeric_limits<T>::min();
            if (isZeroVector) {
                target.zeros();
                return Trv(0);
            }
            target[0] = Trv(2);
            target.template tail<TailLength>(1).zeros();
            return sqrt(sourceNorm0);
        }
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
    __host__ __device__ constexpr void RValueVector<Derived>::static_assert_assign(const Vector auto& target) noexcept {
        using V = std::remove_cvref<decltype(target)>::type;
        constexpr size_t Size1 = SizeAtCompile;
        constexpr size_t Size2 = V::SizeAtCompile;
        static_assert(Size1 == Dynamic || Size2 == Dynamic || Size1 == Size2, "[Error]: Size mismatch between two vector");
        static_assert(V::isComplex || !isComplex, "[Error]: Assign a complex vector to real vector discards imags");
        static_assert(Diffable<V> || !Diffable<This>, "[Error]: Assign a diffable vector to normal vector discards grads");
    }
    /**
     * FIXME: We cannot use unknown references in params of constexpr, refactor once we dump to Clang 20.
     */
    template<class Derived>
    template<Vector V>
    consteval size_t RValueVector<Derived>::maxSizeAtCompile() noexcept {
        return std::max(SizeAtCompile, std::remove_cvref_t<V>::SizeAtCompile);
    }

    template<class Derived>
    template<int GradOrder>
    auto RValueVector<Derived>::grads_impl() const noexcept {
        return GradVector<Derived, GradOrder>(Base::getDerived());
    }

    template<class Derived>
    template<Vector V, ExecutePolicy P>
    void RValueVector<Derived>::assign_for(V& v) const noexcept {
        parallel_for<P>([&, this](size_t i) {
            if constexpr (isReverseDiff)
                v[i] = calc_value(i);
            else
                v[i] = calc(i);
        }, getLength(), 0).wait();
    }

    template<class Derived>
    template<Vector V, ExecutePolicy P, size_t Length>
    void RValueVector<Derived>::assign_simd(V& v) const noexcept {
        using Pack = BestPacket<typename V::ScalarType, Length>::Type;
        constexpr static size_t PacketSize = Pack::size();
        const auto& v0 = Base::getDerived();
        if constexpr (Length != Dynamic) {
            constexpr size_t to = Length / PacketSize * PacketSize;
            for (size_t i = 0; i < to; i += PacketSize)
                v.writePacket(i, v0.template packet<Pack>(i));
            
            constexpr size_t i = Length - Length % PacketSize;
            if constexpr (i != Length) {
                constexpr size_t count = Length - i;
                v.writePacketPartial(i, count, v0.template packetPartial<Pack>(i, count));
            }
        }
        else {
            const size_t length = getLength();
            if (length == 0)
                return;

            const size_t to = length / PacketSize * PacketSize;
            if constexpr (P == Sequential) {
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
                auto future = parallel_for<P>([&, this](size_t i) {
                    const size_t i1 = i * PacketSize;
                    v.writePacket(i1, v0.template packet<Pack>(i1));
                }, numLoop, 0);

                if (to != length) {
                    const size_t count = length % PacketSize;
                    v.writePacketPartial(to, count, v0.template packetPartial<Pack>(to, count));
                }
                future.wait();
            }
        }
    }

    template<class Derived>
    template<Vector V, ExecutePolicy P>
    void RValueVector<Derived>::assign_add_for(V& v) const noexcept {
        parallel_for<P>([&, this](size_t i) {
            if constexpr (isReverseDiff)
                v[i] += calc_value(i);
            else
                v[i] += calc(i);
        }, getLength(), 0).wait();
    }

    template<class Derived>
    template<Vector V, size_t Length>
    void RValueVector<Derived>::assign_add_simd(V& v) const noexcept {
        using Pack = BestPacket<typename V::ScalarType, Length>::Type;
        constexpr static size_t PacketSize = Pack::size();
        const auto& v0 = Base::getDerived();
        if constexpr (Length != Dynamic) {
            constexpr size_t to = Length / PacketSize * PacketSize;
            for (size_t i = 0; i < to; i += PacketSize) {
                const Pack sum = v.template packet<Pack>(i) + v0.template packet<Pack>(i);
                v.writePacket(i, sum);
            }

            constexpr size_t i = Length - Length % PacketSize;
            if constexpr (i != Length) {
                constexpr size_t count = Length - i;
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
}
