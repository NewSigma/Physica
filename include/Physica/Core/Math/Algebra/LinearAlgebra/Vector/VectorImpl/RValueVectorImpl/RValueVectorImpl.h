/*
 * Copyright 2022-2026 Weibo He.
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
#include "Physica/Core/Utils/Range.h"

namespace Physica {
    template<class Derived>
    bool RValueVector<Derived>::operator!=(const Vector auto& other) const noexcept {
        return !(Base::getDerived() == other);
    }

    template<class Derived>
    auto RValueVector<Derived>::operator*(this auto&& self, Scalar auto&& x) noexcept {
        using V = decltype(self);
        using U = decltype(x);
        return VectorExpr<ExprID::Mul, V, U>(std::forward<V>(self), std::forward<U>(x));
    }

    template<class Derived>
    auto RValueVector<Derived>::operator*(this auto&& self, Matrix auto&& m) noexcept {
        using Self = decltype(self);
        using M = decltype(m);
        static_assert(!is_device_obj<M>::value, "[Error]: host-device mismatch");
        return GEVM<Self&&, M&&>(std::forward<Self>(self), std::forward<M>(m));
    }

    template<class Derived>
    auto RValueVector<Derived>::operator-(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (instanceof<Derived, GEMV>)
            return std::forward<Self>(self).getLHS() * (-std::forward<Self>(self).getRHS());
        else
            return VectorExpr<ExprID::Minus, Self>(std::forward<Self>(self));
    }

    template<class Derived>
    template<ExecutePolicy P>
    void RValueVector<Derived>::assign(Vector auto&& v) const noexcept {
        using V = std::remove_cvref<decltype(v)>::type;
        if constexpr (!isDiffable() && v.isDiffable()) {
            Base::getDerived().assign(v.values());
            v.zero_grad();
        }
        else {
            v.assert_assign(Base::getDerived());
            if constexpr (Internal::EnableSIMD<Derived, V>::value && !isReverseDiff()) {
                constexpr static size_t Size = std::max(getSizeAtCompile(), v.getSizeAtCompile());
                assign_simd<V, P, Size>(v);
            }
            else
                assign_for<P>(v);
        }
    }

    template<class Derived>
    template<ExecutePolicy P>
    void RValueVector<Derived>::assign_base(Vector auto&& v) const noexcept {
        assign<P>(v);
    }

    template<class Derived>
    template<ExecutePolicy P>
    void RValueVector<Derived>::assign_add(Vector auto&& v) const noexcept {
        using V = std::remove_cvref<decltype(v)>::type;
        if constexpr (!isDiffable() && v.isDiffable())
            Base::getDerived().assign_add(v.values());
        else {
            v.assert_assign(Base::getDerived());
            if constexpr (Internal::EnableSIMD<Derived, V>::value && !isReverseDiff()) {
                constexpr static size_t Size = std::max(getSizeAtCompile(), v.getSizeAtCompile());
                assign_add_simd<V, Size>(v);
            }
            else
                assign_add_for<P>(v);
        }
    }

    template<class Derived>
    void RValueVector<Derived>::assert_assign(const Vector auto& source) const noexcept {
        static_assert_assign(source);
        if constexpr (std::same_as<Derived, std::remove_cvref_t<decltype(source)>>)
            assert(this != &source && "[Error]: Self assign is likely a bug");
        if constexpr (getSizeAtCompile() == Dynamic && source.getSizeAtCompile() == Dynamic) {
            assert(getLength() == source.getLength() && "[Error]: Size mismatch between two vector");
            assert(!empty());
        }
    }

    template<class Derived>
    void RValueVector<Derived>::assert_assign_lapack(const Vector auto& source) const noexcept {
        static_assert(Internal::EnableLAPACK<Derived, decltype(source)>::value, "[Error]: Invalid expr for LAPACK");
        assert_assign(source);
    }

    template<class Derived>
    decltype(auto) RValueVector<Derived>::calc(size_t index) const noexcept {
        return Base::getDerived().calc(index);
    }

    template<class Derived>
    decltype(auto) RValueVector<Derived>::calc_value(size_t index) const noexcept {
        return Base::getDerived().values().calc(index);
    }

    template<class Derived>
    template<int Size>
    auto RValueVector<Derived>::packet(size_t index) const noexcept -> SIMD<T, Size> {
        using U = typename Derived::ScalarType;
        assert(index + Size <= getLength() && "[Error]: Index out of range");
        if constexpr (Diffable<U>) {
            if constexpr (isForwardDiff()) {
                const auto& x = Base::getDerived();
                auto values = x.values().template packet<Size>(index);
                auto grads = x.grads().template packet<Size>(index);
                return SIMD<T, Size>(std::move(values), std::move(grads));
            }
            else {
                using Uv = U::ValueType;
                Array<Uv, Size> values{};
                for (size_t i = 0; i < Size; ++i, ++index)
                    values[i] = Uv(calc(index));
                SIMD<Uv, Size> packet{};
                packet.load(values.data());
                return SIMD<T, Size>(std::move(packet));
            }
        }
        else {
            Array<U, Size> buffer{};
            for (size_t i = 0; i < Size; ++i, ++index)
                buffer[i] = U(calc(index));
            SIMD<T, Size> packet{};
            packet.load(buffer.data());
            return packet;
        }
    }

    template<class Derived>
    template<int Size>
    auto RValueVector<Derived>::packet(size_t index, size_t count) const noexcept -> SIMD<T, Size> {
        using U = typename Derived::ScalarType;
        assert(index + count <= getLength() && "[Error]: Index out of range");
        if constexpr (Diffable<U>) {
            if constexpr (isForwardDiff()) {
                const auto& x = Base::getDerived();
                auto values = x.values().template packet<Size>(index, count);
                auto grads = x.grads().template packet<Size>(index, count);
                return SIMD<T, Size>(std::move(values), std::move(grads));
            }
            else {
                using Uv = U::ValueType;
                Array<Uv, Size> values{};
                for (size_t i = 0; i < Size; ++i, ++index)
                    values[i] = i < count ? Uv(calc(index)) : Uv(0);
                SIMD<Uv, Size> packet{};
                packet.load(values.data());
                return SIMD<T, Size>(std::move(packet));
            }
        }
        else {
            Array<U, Size> buffer{};
            for (size_t i = 0; i < Size; ++i, ++index)
                buffer[i] = i < count ? U(calc(index)) : U(0);
            SIMD<T, Size> packet{};
            packet.load(buffer.data());
            return packet;
        }
    }

    template<class Derived>
    constexpr auto RValueVector<Derived>::view(this auto&& self) noexcept {
        return View<std::remove_reference_t<decltype(self)>>(self);
    }

    template<class Derived>
    void RValueVector<Derived>::reverse(const Vector auto&, const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff());
        Base::getDerived().reverse(grad);
    }

    template<class Derived>
    void RValueVector<Derived>::resize(const Vector auto& x) {
        resize(x.getLength());
    }

    template<class Derived>
    void RValueVector<Derived>::resize(size_t length) {
        Base::getDerived().resize(length);
    }

    template<class Derived>
    template<size_t Length>
    auto RValueVector<Derived>::head(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        return RVectorBlock<Self, Length>(std::forward<Self>(self), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    auto RValueVector<Derived>::tail(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        return RVectorBlock<Self, Length>(std::forward<Self>(self), from);
    }

    template<class Derived>
    template<size_t Length>
    auto RValueVector<Derived>::segment(this auto&& self, size_t from, size_t to) noexcept {
        using Self = decltype(self);
        return RVectorBlock<Self, Length>(std::forward<Self>(self), from, to);
    }

    template<class Derived>
    template<int Major, size_t Row, size_t Col>
    auto RValueVector<Derived>::reshape(size_t row, size_t col) const noexcept {
        return RValueReshapedVector<const Derived&, Major, Row, Col>(Base::getDerived(), row, col);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto RValueVector<Derived>::reshape_row(size_t row, size_t col) const noexcept {
        return reshape<MatrixMajor::Row, Row, Col>(row, col);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto RValueVector<Derived>::reshape_col(size_t row, size_t col) const noexcept {
        return reshape<MatrixMajor::Col, Row, Col>(row, col);
    }

    template<class Derived>
    auto RValueVector<Derived>::reshape_like(const Matrix auto& mat) const noexcept {
        using M = std::remove_cvref_t<decltype(mat)>;
        constexpr auto Major = MatrixMajor::getMajor<M>();
        if constexpr (Major == MatrixMajor::BothMajor)
            return reshape<MatrixMajor::Col, M::getRowAtCompile(), M::getColAtCompile()>(mat.getRow(), mat.getCol());
        else
            return reshape<Major, M::getRowAtCompile(), M::getColAtCompile()>(mat.getRow(), mat.getCol());
    }

    template<class Derived>
    auto RValueVector<Derived>::transpose(this auto&& self) noexcept {
        using Self = decltype(self);
        return Transpose<Self>(std::forward<Self>(self));
    }

    template<class Derived>
    decltype(auto) RValueVector<Derived>::conjugate(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (isComplex())
            return Conjugate<Self>(std::forward<Self>(self));
        else
            return std::forward<Self>(self);
    }

    template<class Derived>
    auto RValueVector<Derived>::hermite(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (isComplex())
            return Hermite<Self>(std::forward<Self>(self));
        else
            return std::forward<Self>(self).transpose();
    }

    template<class Derived>
    auto RValueVector<Derived>::format() const {
        return FormatedVector<Derived>(Base::getDerived());
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
        return Base::getDerived().squaredNorms().sum();
    }

    template<class Derived>
    auto RValueVector<Derived>::lnSquaredNorm() const -> Tr {
        const auto& derived = Base::getDerived();
        const Tr maxabs = abs(derived).max();
        // We require a small threshold to avoid latter ill-conditioned reciprocal.
        assert(!maxabs.isSubNormal() && "[Error]: Vectors near zero are invalid");
        const Tr factor = reciprocal(maxabs);
        return ln((derived * factor).squaredNorm()) + Tr(2) * ln(maxabs);
    }

    template<class Derived>
    auto RValueVector<Derived>::normInf() const -> Tr {
        return abs(Base::getDerived()).max();
    }

    template<class Derived>
    size_t RValueVector<Derived>::argmax() const noexcept {
        Trv x = std::numeric_limits<Trv>::lowest();
        size_t result = 0;
        for (size_t i = 0; i < getLength(); ++i) {
            Trv y = calc_value(i);
            if (y > x) {
                x = y;
                result = i;
            }
        }
        return result;
    }

    template<class Derived>
    size_t RValueVector<Derived>::argmin() const noexcept {
        Trv x = std::numeric_limits<Trv>::max();
        size_t result = 0;
        for (size_t i = 0; i < getLength(); ++i) {
            Trv y = calc_value(i);
            if (y < x) {
                x = y;
                result = i;
            }
        }
        return result;
    }

    template<class Derived>
    auto RValueVector<Derived>::max() const noexcept -> CoDiff<T> {
        static_assert(!T::isComplex(), "[Error]: Compare between complex number is ill defined");
        assert(getLength() != 0);

        constexpr bool EnableSIMD = Internal::EnableSIMD<Derived>::value;
        if constexpr (isReverseDiff()) {
            size_t index = 0;
            Tv elem = calc_value(0);
            for (size_t i = 1; i < getLength(); ++i) {
                const Tv temp = calc_value(i);
                if (elem < temp) {
                    elem = temp;
                    index = i;
                }
            }
            const auto& result = co_yield calc_value(index);
            calc(index).reverse(result.grad());
        }
        else if constexpr (!EnableSIMD) {
            T result = calc(0);
            for (size_t i = 1; i < getLength(); ++i) {
                T temp = calc(i);
                if (result < temp)
                    result = temp;
            }
            co_return std::move(result);
        }
        else {
            using PacketType = BestPacket<ScalarType, getSizeAtCompile()>::Type;
            constexpr size_t SizeV = getSizeAtCompile();
            constexpr int SizeP = PacketType::size();
            const auto& v = Base::getDerived();
            PacketType buffer(std::numeric_limits<T>::lowest());
            T result;
            if constexpr (SizeV != Dynamic) {
                constexpr size_t to = SizeV / SizeP * SizeP;
                for (size_t i = 0; i < to; i += SizeP)
                    buffer = std::max(v.template packet<SizeP>(i), buffer);
                result = buffer.max();

                constexpr size_t i = SizeV - SizeV % SizeP;
                if constexpr (i != SizeV) {
                    constexpr size_t count = SizeV - i;
                    for (size_t j = 0; j < count; ++j)
                        result = std::max(result, v.calc(i + j));
                }
            }
            else {
                const size_t length = getLength();
                size_t i = 0;
                const size_t to = length / SizeP * SizeP;
                for (; i < to; i += SizeP)
                    buffer = std::max(v.template packet<SizeP>(i), buffer);
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
        static_assert(!T::isComplex(), "[Error]: Compare between complex number is ill defined");
        assert(getLength() != 0);

        constexpr bool EnableSIMD = Internal::EnableSIMD<Derived>::value;
        if constexpr (isReverseDiff()) {
            size_t index = 0;
            Tv elem = calc_value(0);
            for (size_t i = 1; i < getLength(); ++i) {
                const Tv temp = calc_value(i);
                if (elem > temp) {
                    elem = temp;
                    index = i;
                }
            }
            const auto& result = co_yield calc_value(index);
            calc(index).reverse(result.grad());
        }
        else if constexpr (EnableSIMD) {
            using PacketType = BestPacket<ScalarType, getSizeAtCompile()>::Type;
            constexpr size_t SizeV = getSizeAtCompile();
            constexpr int SizeP = PacketType::size();
            const auto& v = Base::getDerived();
            PacketType buffer(std::numeric_limits<T>::max());
            T result;
            if constexpr (SizeV != Dynamic) {
                constexpr size_t to = SizeV / SizeP * SizeP;
                for (size_t i = 0; i < to; i += SizeP)
                    buffer = std::min(v.template packet<SizeP>(i), buffer);
                result = buffer.min();

                constexpr size_t i = SizeV - SizeV % SizeP;
                if constexpr (i != SizeV) {
                    constexpr size_t count = SizeV - i;
                    for (size_t j = 0; j < count; ++j)
                        result = std::min(result, v.calc(i + j));
                }
            }
            else {
                const size_t length = getLength();
                size_t i = 0;
                const size_t to = length / SizeP * SizeP;
                for (; i < to; i += SizeP)
                    buffer = std::min(v.template packet<SizeP>(i), buffer);
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
            for (size_t i = 1; i < getLength(); ++i) {
                T temp = calc(i);
                if (result > temp)
                    result = temp;
            }
            co_return std::move(result);
        }
    }

    template<class Derived>
    auto RValueVector<Derived>::sum(this auto&& self) noexcept -> CoDiff<T> {
        assert(self.getLength() != 0 && "[Error]: Sum of a empty vector is not well defined");
        if constexpr (isReverseDiff()) {
            using Self = decltype(self);
            decltype(auto) v = decay_rvalue(std::forward<Self>(self));
            auto& result = co_yield v.values().sum();
            v.reverse(result.grad());
        }
        else if constexpr (Internal::EnableSIMD<Derived>::value) {
            using PacketType = BestPacket<ScalarType, getSizeAtCompile()>::Type;
            constexpr size_t SizeV = getSizeAtCompile();
            constexpr int SizeP = PacketType::size();
            auto buffer = PacketType::zeros();
            if constexpr (SizeV != Dynamic) {
                constexpr size_t to = SizeV / SizeP * SizeP;
                for (size_t i = 0; i < to; i += SizeP)
                    buffer += self.template packet<SizeP>(i);

                constexpr size_t i = SizeV - SizeV % SizeP;
                if constexpr (i != SizeV) {
                    constexpr size_t count = SizeV - i;
                    buffer += self.template packet<SizeP>(i, count);
                }
            }
            else {
                const size_t length = self.getLength();
                size_t i = 0;
                const size_t to = length / SizeP * SizeP;
                for (; i < to; i += SizeP)
                    buffer += self.template packet<SizeP>(i);

                if (to != length) {
                    const size_t count = length - i;
                    buffer += self.template packet<SizeP>(i, count);
                }
            }
            co_return buffer.sum();
        }
        else {
            auto result = T(0);
            for (size_t i = 0; i < self.getLength(); ++i)
                result += self.calc(i);
            co_return std::move(result);
        }
    }

    template<class Derived>
    auto RValueVector<Derived>::mean() const noexcept -> CoDiff<T> {
        return Base::getDerived().sum() / Trv(getLength());
    }

    template<class Derived>
    auto RValueVector<Derived>::mean_stable() const noexcept -> CoDiff<T> {
        auto result = T(0);
        const auto& v = Base::getDerived();
        for (size_t i = 0; i < getLength(); ++i)
            result.toNextMean(i, v.calc(i));

        if constexpr (isReverseDiff()) {
            auto& y = co_yield std::move(result);
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
        if constexpr (isReverseDiff()) {
            auto& tmp = co_yield result.value();
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
    auto RValueVector<Derived>::lnSumExp(this const auto& self) noexcept -> CoDiff<T> {
        if constexpr (self.isLValueVector()) {
            Tv m;
            if constexpr (isComplex())
                m = self.values().reals().max();
            else
                m = self.values().max();

            auto y = ln(exp(self - m).sum() + Trv(std::numeric_limits<Trv>::min())) + m;
            if constexpr (isReverseDiff()) {
                auto& tmp = co_yield y.value();
                y.reverse_final(tmp.grad());
            }
            else
                co_return std::move(y);
        }
        else
            co_return DenseVector<T, getSizeAtCompile(), HostAllocator<T>>(self).lnSumExp();
    }

    template<class Derived>
    auto RValueVector<Derived>::crossEntropy(size_t index) const noexcept -> CoDiff<T> {
        assert(index < getLength() && "[Error]: Index overflow");
        auto y = Base::getDerived().lnSumExp() - calc_value(index);
        if constexpr (isReverseDiff()) {
            auto& tmp = co_yield y.value();
            auto g = tmp.grad();
            y.reverse(g);
            calc(index).reverse(-g);
        }
        else
            co_return std::move(y);
    }

    template<class Derived>
    auto RValueVector<Derived>::lnSoftmax(size_t index) const noexcept -> CoDiff<T> {
        assert(index < getLength() && "[Error]: Index overflow");
        auto y = -crossEntropy(index);
        if constexpr (isReverseDiff()) {
            auto& tmp = co_yield y.value();
            y.reverse(tmp.grad());
        }
        else
            co_return std::move(y);
    }

    template<class Derived>
    auto RValueVector<Derived>::softmax(size_t index) const noexcept -> CoDiff<T> {
        assert(index < getLength() && "[Error]: Index overflow");
        auto y = exp(lnSoftmax(index));
        if constexpr (isReverseDiff()) {
            auto& tmp = co_yield y.value();
            y.reverse_final(tmp.grad());
        }
        else
            co_return std::move(y);
    }

    template<class Derived>
    auto RValueVector<Derived>::prod() const noexcept -> CoDiff<T> {
        assert(getLength() != 0);
        if constexpr (isReverseDiff()) {
            Tv result = calc_value(0);
            for (size_t i = 1; i < getLength(); ++i)
                result *= calc_value(i);

            auto& tmp = co_yield std::move(result);
            const auto& v = Base::getDerived();
            v.reverse(reciprocal(v.values()) * (tmp.value() * tmp.grad()));
        }
        else {
            auto result = calc(0);
            for (size_t i = 1; i < getLength(); ++i)
                result *= calc(i);
            co_return std::move(result);
        }
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
    auto RValueVector<Derived>::householder(Vector auto& __restrict target) const __restrict -> Tr {
        constexpr size_t Size = std::max(getSizeAtCompile(), target.getSizeAtCompile());
        constexpr size_t TailLength = Size > 0 ? (Size - 1) : Dynamic;
        assert(getLength() == target.getLength());
        assert(getLength() > 1 && "[Error]: Unnecessary householder call");

        const T v0 = calc(0);
        const Tr sourceNorm0 = v0.squaredNorm();
        const Tr squaredTailNorm = Base::getDerived().template tail<TailLength>(1).squaredNorm();
        if (!squaredTailNorm.isSubNormal()) [[likely]] {
            const Tr norm = sqrt(squaredTailNorm + sourceNorm0);
            target[0] = Tr(1) + abs(v0) / norm;
            target.template tail<TailLength>(1) = Base::getDerived().template tail<TailLength>(1) * reciprocal(v0 + unit(v0.value()) * norm);
            return norm;
        }

        const bool isZeroVector = sourceNorm0.isSubNormal();
        if (isZeroVector) {
            target.zeros();
            return Trv(0);
        }
        target[0] = Trv(2);
        target.template tail<TailLength>(1).zeros();
        return sqrt(sourceNorm0);
    }

    template<class Derived>
    decltype(auto) RValueVector<Derived>::reals(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (isComplex())
            return RealVectorR<Self>(std::forward<Self>(self));
        else
            return std::forward<Self>(self);
    }

    template<class Derived>
    auto RValueVector<Derived>::imags(this auto&& self) noexcept {
        using Self = decltype(self);
        return ImagVectorR<Self>(std::forward<Self>(self));
    }

    template<class Derived>
    auto RValueVector<Derived>::squaredNorms(this auto&& self) noexcept {
        using Self = decltype(self);
        return SquaredNormVector<Self>(std::forward<Self>(self));
    }

    template<class Derived>
    auto RValueVector<Derived>::norms(this auto&& self) noexcept {
        using Self = decltype(self);
        return NormVector<Self>(std::forward<Self>(self));
    }

    template<class Derived>
    decltype(auto) RValueVector<Derived>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (isDiffable())
            return ValueVector<Self>(std::forward<Self>(self));
        else
            return std::forward<Self>(self);
    }

    template<class Derived>
    template<int GradOrder>
    auto RValueVector<Derived>::grads(this auto&& self) noexcept {
        using Self = decltype(self);
        return GradVector<Self, GradOrder>(std::forward<Self>(self));
    }

    template<class Derived>
    template<int MaskOrder>
    decltype(auto) RValueVector<Derived>::grads_mask(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (MaskOrder < T::Order)
            return GradMaskVector<Self, MaskOrder>(std::forward<Self>(self));
        else
            return std::forward<Self>(self);
    }

    template<class Derived>
    bool RValueVector<Derived>::isZero() const {
        for (size_t i = 0; i < getLength(); ++i)
            if (!calc(i).isZero())
                return false;
        return true;
    }

    template<class Derived>
    bool RValueVector<Derived>::isFinite() const {
        for (size_t i = 0; i < getLength(); ++i)
            if (!calc(i).isFinite())
                return false;
        return true;
    }

    template<class Derived>
    bool RValueVector<Derived>::isSubNormal() const {
        for (size_t i = 0; i < getLength(); ++i)
            if (!calc(i).isSubNormal())
                return false;
        return true;
    }

    template<class Derived>
    __host__ __device__ consteval bool RValueVector<Derived>::isComplex() noexcept {
        return ScalarType::isComplex();
    }

    template<class Derived>
    __host__ __device__ consteval bool RValueVector<Derived>::isDiffable() noexcept {
        return ScalarType::isDiffable();
    }

    template<class Derived>
    __host__ __device__ consteval bool RValueVector<Derived>::isForwardDiff() noexcept {
        return ScalarType::isForwardDiff();
    }

    template<class Derived>
    __host__ __device__ consteval bool RValueVector<Derived>::isReverseDiff() noexcept {
        return ScalarType::isReverseDiff();
    }

    template<class Derived>
    __host__ __device__ consteval bool RValueVector<Derived>::isLValueVector() noexcept {
        return false;
    }

    template<class Derived>
    __host__ __device__ consteval bool RValueVector<Derived>::isCompact() noexcept {
        return requires{ std::declval<Derived>().data(); };
    }

    template<class Derived>
    __host__ __device__ consteval bool RValueVector<Derived>::isSparse() noexcept {
        return requires{ std::declval<Derived>().getNumNonzero(); };
    }

    template<class Derived>
    __host__ __device__ consteval bool RValueVector<Derived>::isFastAssign() noexcept {
        return false;
    }

    template<class Derived>
    __host__ __device__ consteval bool RValueVector<Derived>::isFastPacket() noexcept {
        return false;
    }

    template<class Derived>
    __host__ __device__ consteval auto RValueVector<Derived>::getSizeAtCompile() noexcept {
        return Derived::getSizeAtCompile();
    }

    template<class Derived>
    __host__ __device__ consteval void RValueVector<Derived>::static_assert_assign(const Scalar auto& source) noexcept {
        using U = std::remove_cvref<decltype(source)>::type;
        T::template static_assert_assign<U>();
    }

    template<class Derived>
    __host__ __device__ consteval void RValueVector<Derived>::static_assert_assign(const Vector auto& source) noexcept {
        using Src = std::remove_cvref<decltype(source)>::type;
        constexpr size_t Size1 = getSizeAtCompile();
        constexpr size_t Size2 = source.getSizeAtCompile();
        static_assert(Size1 == Dynamic || Size2 == Dynamic || Size1 == Size2, "[Error]: Size mismatch between two vector");

        using U = Src::ScalarType;
        T::template static_assert_assign<U>();
    }

    template<class Derived>
    consteval int RValueVector<Derived>::calcBlockingSize(int CacheSize) noexcept {
        int result = 1;
        while (result * int(sizeof(Trv)) < CacheSize)
            result *= 2;
        return result / 2;
    }

    template<class Derived>
    template<ExecutePolicy P>
    void RValueVector<Derived>::assign_for(Vector auto& v) const noexcept {
        auto it = zip(v.view(), Base::getDerived().view()).begin();
        parallel_for<P>([it](size_t i) noexcept {
            auto [lhs, rhs] = it + i;
            *lhs = *rhs;
        }, getLength(), 0).wait();
    }

    template<class Derived>
    template<Vector V, ExecutePolicy P, size_t Length>
    void RValueVector<Derived>::assign_simd(V& v) const noexcept {
        constexpr int Size = BestPacket<typename V::ScalarType, Length>::Size;
        auto it = zip(v.view(), Base::getDerived().view()).begin();
        if constexpr (Length != Dynamic) {
            constexpr size_t to = Length / Size * Size;
            for (size_t i = 0; i < to; i += Size) {
                auto [lhs, rhs] = it + i;
                lhs.store(rhs.template load<Size>());
            }

            for (size_t i = Length - Length % Size; i < Length; ++i) {
                auto [lhs, rhs] = it + i;
                *lhs = *rhs;
            }
        }
        else {
            const size_t length = getLength();
            const size_t to = length / Size * Size;
            if constexpr (P == Sequential) {
                size_t i = 0;
                for (; i < to; i += Size) {
                    auto [lhs, rhs] = it + i;
                    lhs.store(rhs.template load<Size>());
                }

                for (; i < length; ++i) {
                    auto [lhs, rhs] = it + i;
                    *lhs = *rhs;
                }
            }
            else {
                auto future = parallel_for<P>([it](size_t i) {
                    auto [lhs, rhs] = it + i * Size;
                    lhs.store(rhs.template load<Size>());
                }, to / Size, 0);

                for (size_t i = to; i < length; ++i) {
                    auto [lhs, rhs] = it + i;
                    *lhs = *rhs;
                }
                future.wait();
            }
        }
    }

    template<class Derived>
    template<ExecutePolicy P>
    void RValueVector<Derived>::assign_add_for(Vector auto& v) const noexcept {
        auto it = zip(v.view(), Base::getDerived().view()).begin();
        parallel_for<P>([it](size_t i) noexcept {
            auto [lhs, rhs] = it + i;
            *lhs += *rhs;
        }, getLength(), 0).wait();
    }

    template<class Derived>
    template<Vector V, size_t Length>
    void RValueVector<Derived>::assign_add_simd(V& v) const noexcept {
        constexpr int Size = BestPacket<typename V::ScalarType, Length>::Size;
        auto it = zip(v.view(), Base::getDerived().view()).begin();
        if constexpr (Length != Dynamic) {
            constexpr size_t to = Length / Size * Size;
            for (size_t i = 0; i < to; i += Size) {
                auto [lhs, rhs] = it + i;
                lhs.store(lhs.template load<Size>() + rhs.template load<Size>());
            }

            for (size_t i = Length - Length % Size; i < Length; ++i) {
                auto [lhs, rhs] = it + i;
                *lhs += *rhs;
            }
        }
        else {
            const size_t length = v.getLength();
            const size_t to = length / Size * Size;
            size_t i = 0;
            for (; i < to; i += Size) {
                auto [lhs, rhs] = it + i;
                lhs.store(lhs.template load<Size>() + rhs.template load<Size>());
            }

            for (; i < length; ++i) {
                auto [lhs, rhs] = it + i;
                *lhs += *rhs;
            }
        }
    }
}
