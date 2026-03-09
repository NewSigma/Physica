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

#include "../LValueVector.h"

namespace Physica {
    template<class Derived>
    auto LValueVector<Derived>::operator=(Scalar auto x) noexcept -> Derived& {
        Base::static_assert_assign(x);
        if constexpr (!std::same_as<T, std::remove_cvref_t<decltype(x)>>)
            return operator=(T(x));
        else {
            for (size_t i = 0; i < Base::getLength(); ++i)
                (*this)[i] = x;
            return Base::getDerived();
        }
    }

    template<class Derived>
    void LValueVector<Derived>::operator+=(Scalar auto x) noexcept {
        auto& v = Base::getDerived();
        (v + x).assign(v);
    }

    template<class Derived>
    void LValueVector<Derived>::operator-=(Scalar auto x) noexcept {
        auto& v = Base::getDerived();
        (v - x).assign(v);
    }

    template<class Derived>
    void LValueVector<Derived>::operator*=(Scalar auto x) noexcept {
        auto& v = Base::getDerived();
        (v * x).assign(v);
    }

    template<class Derived>
    void LValueVector<Derived>::operator/=(Scalar auto x) noexcept {
        auto& v = Base::getDerived();
        (v / x).assign(v);
    }

    template<class Derived>
    template<ExecutePolicy P>
    Derived& LValueVector<Derived>::operator=(const Vector auto& v) {
        if constexpr (std::is_same<const Derived&, decltype(v)>::value)
            assert(this != &v && "[Error]: Self assign is likely a bug");
        Derived& x = Base::getDerived();
        x.resize(v);
        x.assert_assign(v);
        v.template assign<P>(x);
        return x;
    }

    template<class Derived>
    void LValueVector<Derived>::operator+=(const Vector auto& v) {
        v.assign_add(Base::getDerived());
    }

    template<class Derived>
    void LValueVector<Derived>::operator-=(const Vector auto& v) {
        Base::getDerived() += -v;
    }

    template<class Derived>
    decltype(auto) LValueVector<Derived>::operator[](this auto&& self, size_t index) {
        return *self.data_ptr(index);
    }

    template<class Derived>
    decltype(auto) LValueVector<Derived>::calc(size_t index) const {
        return operator[](index);
    }

    template<class Derived>
    auto LValueVector<Derived>::calc_value(size_t index) const -> Tv {
        return operator[](index).value();
    }

    template<class Derived>
    void LValueVector<Derived>::writePacket(const Packet auto packet, size_t index) noexcept {
        using Pack = std::remove_cvref_t<decltype(packet)>;
        using U = Traits<Pack>::ScalarType;
        if constexpr (U::isForwardDiff) {
            for (size_t i = 0; i < Pack::size(); ++i, ++index)
                (*this)[index] = packet[i];
        }
        else {
            Array<U, Pack::size()> buffer{};
            packet.store(buffer.data());
            for (size_t i = 0; i < Pack::size(); ++i, ++index)
                (*this)[index] = buffer[i];
        }
    }

    template<class Derived>
    void LValueVector<Derived>::writePacket(const Packet auto packet, size_t index, size_t count) noexcept {
        using Pack = std::remove_cvref_t<decltype(packet)>;
        using U = Traits<Pack>::ScalarType;
        if constexpr (U::isForwardDiff) {
            for (size_t i = 0; i < count; ++i, ++index)
                (*this)[index] = packet[i];
        }
        else {
            Array<U, Pack::size()> buffer{};
            packet.store(buffer.data());
            for (size_t i = 0; i < count; ++i, ++index)
                (*this)[index] = buffer[i];
        }
    }

    template<class Derived>
    constexpr auto LValueVector<Derived>::view(this auto&& self) noexcept {
        return LVectorView<std::remove_reference_t<decltype(self)>>(self);
    }

    template<class Derived>
    auto LValueVector<Derived>::sum() const -> CoDiff<T> {
        if constexpr (isReverseDiff) {
            auto& result = co_yield Base::getDerived().values().sum();
            const auto& grad = result.grad();
            for (size_t i = 0; i < Base::getLength(); ++i)
                (*this)[i].reverse(grad);
        }
        else
            co_return Base::sum();
    }

    template<class Derived>
    void LValueVector<Derived>::toNextMean(size_t lastNumSample, const Vector auto& sample) noexcept {
        Base::assert_assign(sample);
        auto& mean = Base::getDerived();
        if (lastNumSample == 0) [[unlikely]] {
            sample.assign(mean);
            return;
        }

        const T factor1 = T(lastNumSample);
        const T factor2 = reciprocal(T(lastNumSample + 1));
        mean = (factor1 * mean + sample) * factor2;
    }

    template<class Derived>
    void LValueVector<Derived>::toNextVariance(Derived& mean, size_t lastNumSample, const Vector auto& sample) noexcept {
        Base::assert_assign(sample);
        auto& var = Base::getDerived();
        if (lastNumSample == 0) {
            sample.assign(mean);
            var.zeros();
            return;
        }

        const T factor1 = T(lastNumSample);
        const T factor2 = reciprocal(T(lastNumSample + 1));
        var = (var + square(mean - sample) * factor2) * (factor1 * factor2);
        mean.toNextMean(lastNumSample, sample);
    }

    template<class Derived>
    void LValueVector<Derived>::reverse(const auto& grad) const noexcept {
        using U = std::remove_cvref_t<decltype(grad)>;
        static_assert(std::same_as<typename T::GradType, typename U::ScalarType>, "[Error]: Inconsistent ScalarType");
        static_assert(isReverseDiff);
        if constexpr (Scalar<U>)
            Base::getConstCastDerived().grads() += grad;
        else if constexpr (Vector<U>) {
            assert(Base::getLength() == grad.getLength());
            grad.assign_add(Base::getConstCastDerived().grads());
        }
        else {
            static_assert(Matrix<U>, "[Error]: Unexpected type");
            assert(Base::getLength() == grad.getRow());
            reverse(grad.sum_cols());
        }
    }

    template<class Derived>
    template<size_t Length>
    auto LValueVector<Derived>::head(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        return LVectorBlock<Self, Length>(std::forward<Self>(self), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    auto LValueVector<Derived>::tail(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        return LVectorBlock<Self, Length>(std::forward<Self>(self), from);
    }

    template<class Derived>
    template<size_t Length>
    auto LValueVector<Derived>::segment(this auto&& self, size_t from, size_t to) noexcept {
        using Self = decltype(self);
        return LVectorBlock<Self, Length>(std::forward<Self>(self), from, to);
    }

    template<class Derived>
    void LValueVector<Derived>::zeros() noexcept {
        for (size_t i = 0; i < Base::getLength(); ++i)
            (*this)[i] = Trv(0);
    }

    template<class Derived>
    void LValueVector<Derived>::clamp_min(Tv minimum) {
        for (size_t i = 0; i < Base::getLength(); ++i)
            (*this)[i] = std::max((*this)[i], minimum);
    }

    template<class Derived>
    void LValueVector<Derived>::clamp_max(Tv maximum) {
        for (size_t i = 0; i < Base::getLength(); ++i)
            (*this)[i] = std::min((*this)[i], maximum);
    }

    template<class Derived>
    void LValueVector<Derived>::toUnit() {
        auto& x = Base::getDerived();
        x *= reciprocal(x.norm());
    }

    template<class Derived>
    void LValueVector<Derived>::standardize() {
        assert(Base::getLength() > 1);
        auto& x = Base::getDerived();
        const T mean = Base::mean();
        const T factor = reciprocal(Base::deviation());
        x = (x - mean) * factor;
    }

    template<class Derived>
    auto LValueVector<Derived>::householder() -> Tr {
        return householder(Base::getDerived());
    }

    template<class Derived>
    template<RNG R>
    void LValueVector<Derived>::random_uniform() {
        for (size_t i = 0; i < this->getLength(); ++i)
            this->operator[](i) = T::template random_uniform<R>();
    }

    template<class Derived>
    template<RNG R>
    void LValueVector<Derived>::random_normal() {
        for (size_t i = 0; i < this->getLength(); ++i)
            this->operator[](i) = T::template random_normal<R>();
    }

    template<class Derived>
    template<RNG R>
    void LValueVector<Derived>::random_any(auto& distribution) {
        for (size_t i = 0; i < this->getLength(); ++i)
            this->operator[](i) = T::template random_any<R>(distribution);
    }

    template<class Derived>
    decltype(auto) LValueVector<Derived>::reals() noexcept {
        if constexpr (isComplex)
            return RealVectorL<Derived>(Base::getDerived());
        else
            return Base::getDerived();
    }

    template<class Derived>
    decltype(auto) LValueVector<Derived>::reals() const noexcept {
        return Base::getConstCastDerived().reals();
    }

    template<class Derived>
    template<int GradOrder>
    auto LValueVector<Derived>::grads() const noexcept {
        return Base::template grads_impl<GradOrder>();
    }
    /**
     * Add this function because we cannot simply return &(*this)[index], it is invalid to dereference a device pointer on host.
     */
    template<class Derived>
    auto LValueVector<Derived>::data_ptr(this auto&& self, size_t index) noexcept {
        assert(index < self.getLength() && "[Error]: Index out of range");
        return self.getDerived().data_ptr(index);
    }
}
