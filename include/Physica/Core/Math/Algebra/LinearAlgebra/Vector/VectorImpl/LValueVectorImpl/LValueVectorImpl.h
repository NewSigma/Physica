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

#include "../LValueVector.h"

namespace Physica {
    template<class Derived>
    LValueVector<Derived>& LValueVector<Derived>::operator=(const This& v) {
        return operator=<Derived>(v);
    }

    template<class Derived>
    LValueVector<Derived>& LValueVector<Derived>::operator=(This&& v) {
        return operator=<Derived>(v);
    }

    template<class Derived>
    template<ExecutePolicy P>
    Derived& LValueVector<Derived>::operator=(const Vector auto& v) {
        if constexpr (std::is_same<const Derived&, decltype(v)>::value)
            assert(this != &v && "[Error]: Self assign is likely a bug");
        Derived& x = Base::getDerived();
        x.resize(v.getLength());
        v.template assign<P>(x);
        return x;
    }

    template<class Derived>
    Derived& LValueVector<Derived>::operator=(const Scalar auto& x) {
        constexpr bool isReverseDiffX = ReverseDiff<decltype(x)>;
        static_assert(!isReverseDiff || !isReverseDiffX, "[Error]: Assign a diffable scalar to diffable vector discards grads");
        for (size_t i = 0; i < Base::getLength(); ++i) {
            if constexpr (isReverseDiffX)
                (*this)[i] = x.value();
            else
                (*this)[i] = x;
        }
        return Base::getDerived();
    }

    template<class Derived>
    void LValueVector<Derived>::operator+=(const Vector auto& v) {
        v.assign_add(Base::getDerived());
    }

    template<class Derived>
    void LValueVector<Derived>::operator-=(const Vector auto& v) {
        Base::getDerived() += -v; // To avoid alias
    }

    template<class Derived>
    auto LValueVector<Derived>::operator[](size_t index) -> RefTy {
        return *data_ptr(index);
    }

    template<class Derived>
    auto LValueVector<Derived>::operator[](size_t index) const -> ConstRefTy {
        return *data_ptr(index);
    }

    template<class Derived>
    auto LValueVector<Derived>::calc(size_t index) const -> ConstRefTy {
        return operator[](index);
    }

    template<class Derived>
    auto LValueVector<Derived>::calc_value(size_t index) const -> Tv {
        return operator[](index).value();
    }

    template<class Derived>
    void LValueVector<Derived>::writePacket(size_t index, const Packet auto packet) {
        using Pack = std::remove_cvref_t<decltype(packet)>;
        using T = Traits<Pack>::ScalarType;
        if constexpr (T::isForwardDiff) {
            for (size_t i = 0; i < Pack::size(); ++i, ++index)
                (*this)[index] = packet[i];
        }
        else {
            ScalarType buffer[Pack::size()];
            packet.store(buffer);
            for (size_t i = 0; i < Pack::size(); ++i, ++index)
                (*this)[index] = buffer[i];
        }
    }

    template<class Derived>
    void LValueVector<Derived>::writePacketPartial(size_t index, size_t count, const Packet auto packet) {
        using Pack = std::remove_cvref_t<decltype(packet)>;
        using T = Traits<Pack>::ScalarType;
        if constexpr (T::isForwardDiff) {
            for (size_t i = 0; i < count; ++i, ++index)
                (*this)[index] = packet[i];
        }
        else {
            ScalarType buffer[Pack::size()];
            packet.store(buffer);
            for (size_t i = 0; i < count; ++i, ++index)
                (*this)[index] = buffer[i];
        }
    }

    template<class Derived>
    auto LValueVector<Derived>::sum() const -> CoDiff<ScalarType> {
        if constexpr (isReverseDiff) {
            auto result = co_yield Base::values().sum();
            const auto& grad = result.grad();
            for (size_t i = 0; i < Base::getLength(); ++i)
                (*this)[i].reverse(grad);
        }
        else
            co_return Base::sum();
    }

    template<class Derived>
    void LValueVector<Derived>::reverse(const auto& grad) const noexcept requires(isReverseDiff) {
        using U = std::remove_cvref_t<decltype(grad)>;
        static_assert(std::same_as<typename ScalarType::GradType, typename U::ScalarType>, "[Error]: Inconsistent ScalarType");
        if constexpr (Scalar<U>) {
            for (size_t i = 0; i < Base::getLength(); ++i)
                (*this)[i].reverse(grad);
        }
        else {
            static_assert(Vector<U>, "[Error]: Unexpected type");
            assert(Base::getLength() == grad.getLength());
            for (size_t i = 0; i < Base::getLength(); ++i)
                (*this)[i].reverse(grad.calc(i));
        }
    }

    template<class Derived>
    template<size_t Length>
    auto LValueVector<Derived>::head(size_t to) noexcept {
        return BlockType<Length>(Base::getDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    const auto LValueVector<Derived>::head(size_t to) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    auto LValueVector<Derived>::tail(size_t from) noexcept {
        return BlockType<Length>(Base::getDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    const auto LValueVector<Derived>::tail(size_t from) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    auto LValueVector<Derived>::segment(size_t from, size_t to) noexcept {
        return BlockType<Length>(Base::getDerived(), from, to);
    }

    template<class Derived>
    template<size_t Length>
    const auto LValueVector<Derived>::segment(size_t from, size_t to) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), from, to);
    }

    template<class Derived>
    void LValueVector<Derived>::clamp_min(const Tv& minimum) {
        for (size_t i = 0; i < Base::getLength(); ++i)
            (*this)[i] = std::max((*this)[i], minimum);
    }

    template<class Derived>
    void LValueVector<Derived>::clamp_max(const Tv& maximum) {
        for (size_t i = 0; i < Base::getLength(); ++i)
            (*this)[i] = std::min((*this)[i], maximum);
    }

    template<class Derived>
    void LValueVector<Derived>::toUnit() {
        auto& x = Base::getDerived();
        x *= reciprocal(x.norm());
    }

    template<class Derived>
    void LValueVector<Derived>::normalize() {
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
            this->operator[](i) = ScalarType::template random_uniform<R>();
    }

    template<class Derived>
    template<RNG R>
    void LValueVector<Derived>::random_normal() {
        for (size_t i = 0; i < this->getLength(); ++i)
            this->operator[](i) = ScalarType::template random_normal<R>();
    }

    template<class Derived>
    template<RNG R>
    void LValueVector<Derived>::random_any(auto& distribution) {
        for (size_t i = 0; i < this->getLength(); ++i)
            this->operator[](i) = ScalarType::template random_any<R>(distribution);
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
    auto LValueVector<Derived>::data_ptr(size_t index) noexcept -> PtrTy {
        assert(index < Base::getLength() && "[Error]: Index out of range");
        return Base::getDerived().data_ptr(index);
    }

    template<class Derived>
    auto LValueVector<Derived>::data_ptr(size_t index) const noexcept -> ConstPtrTy {
        return const_cast<This&>(*this).data_ptr(index);
    }
}
