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

#include "../LValueVector.h"

namespace Physica::Core {
    template<class Derived>
    inline LValueVector<Derived>& LValueVector<Derived>::operator=(const This& v) {
        return operator=<Derived>(v);
    }

    template<class Derived>
    inline LValueVector<Derived>& LValueVector<Derived>::operator=(This&& v) {
        return operator=<Derived>(v);
    }

    template<class Derived>
    template<Vector V, class Executor>
    inline Derived& LValueVector<Derived>::operator=(const V& v) {
        if constexpr (std::is_same<Derived, V>::value)
            assert(this != &v && "[Error]: Self assign is likely a bug");
        Derived& v1 = Base::getDerived();
        v1.resize(v.getLength());
        v.template assignTo<Derived, Executor>(v1);
        return v1;
    }

    template<class Derived>
    template<Scalar T>
    inline Derived& LValueVector<Derived>::operator=(const T& x) {
        const ScalarType x1 = x;
        for (size_t i = 0; i < Base::getLength(); ++i)
            (*this)[i] = x1;
        return Base::getDerived();
    }

    template<class Derived>
    inline LValueVector<Derived>::RefTy LValueVector<Derived>::operator[](size_t index) {
        return *data_ptr(index);
    }

    template<class Derived>
    inline LValueVector<Derived>::ConstRefTy LValueVector<Derived>::operator[](size_t index) const {
        return const_cast<This&>(*this).operator[](index);
    }

    template<class Derived>
    LValueVector<Derived>::ScalarType LValueVector<Derived>::calc(size_t index) const {
        return ScalarType(*data_ptr(index));
    }

    template<class Derived>
    template<class AnyPacket>
    void LValueVector<Derived>::writePacket(size_t index, const AnyPacket packet) {
        using T = Traits<AnyPacket>::ScalarType;
        if constexpr (T::isForwardDiff) {
            for (size_t i = 0; i < AnyPacket::size(); ++i, ++index)
                (*this)[index] = packet[i];
        }
        else {
            ScalarType buffer[AnyPacket::size()];
            packet.store(buffer);
            for (size_t i = 0; i < AnyPacket::size(); ++i, ++index)
                (*this)[index] = buffer[i];
        }
    }

    template<class Derived>
    template<class AnyPacket>
    void LValueVector<Derived>::writePacketPartial(size_t index, size_t count, const AnyPacket packet) {
        using T = Traits<AnyPacket>::ScalarType;
        if constexpr (T::isForwardDiff) {
            for (size_t i = 0; i < count; ++i, ++index)
                (*this)[index] = packet[i];
        }
        else {
            ScalarType buffer[AnyPacket::size()];
            packet.store(buffer);
            for (size_t i = 0; i < count; ++i, ++index)
                (*this)[index] = buffer[i];
        }
    }

    template<class Derived>
    CoDiff<typename LValueVector<Derived>::ScalarType> LValueVector<Derived>::sum() const {
        if constexpr (isReverseDiff) {
            auto result = co_yield toValueVector(Base::getDerived()).sum();
            const auto& grad = result.grad();
            for (size_t i = 0; i < Base::getLength(); ++i)
                (*this)[i].reverse(grad);
        }
        else
            co_return Base::sum();
    }

    template<class Derived>
    template<Vector T>
    void LValueVector<Derived>::reverse(const T& grad) const noexcept requires(isReverseDiff) {
        using GradType = ScalarType::GradType;
        static_assert(std::is_same<GradType, typename T::ScalarType>::value, "[Error]: Inconsistent ScalarType");
        for (size_t i = 0; i < Base::getLength(); ++i)
            (*this)[i].reverse(grad.calc(i));
    }

    template<class Derived>
    template<size_t Length>
    inline auto LValueVector<Derived>::head(size_t to) noexcept {
        return BlockType<Length>(Base::getDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    inline const auto LValueVector<Derived>::head(size_t to) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    inline auto LValueVector<Derived>::tail(size_t from) noexcept {
        return BlockType<Length>(Base::getDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    inline const auto LValueVector<Derived>::tail(size_t from) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    inline auto LValueVector<Derived>::segment(size_t from, size_t to) noexcept {
        return BlockType<Length>(Base::getDerived(), from, to);
    }

    template<class Derived>
    template<size_t Length>
    inline const auto LValueVector<Derived>::segment(size_t from, size_t to) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), from, to);
    }

    template<class Derived>
    inline void LValueVector<Derived>::toUnit() {
        Base::getDerived() *= reciprocal(Base::getDerived().norm());
    }

    template<class Derived>
    template<RandomGenerator R>
    inline void LValueVector<Derived>::random_uniform() {
        for (size_t i = 0; i < this->getLength(); ++i)
            this->operator[](i) = ScalarType::template random_uniform<R>();
    }

    template<class Derived>
    template<RandomGenerator R>
    inline void LValueVector<Derived>::random_normal() {
        for (size_t i = 0; i < this->getLength(); ++i)
            this->operator[](i) = ScalarType::template random_normal<R>();
    }

    template<class Derived>
    template<class Distribution, RandomGenerator R>
    inline void LValueVector<Derived>::random_any(Distribution& dist) {
        for (size_t i = 0; i < this->getLength(); ++i)
            this->operator[](i) = ScalarType::template random_any<Distribution, R>(dist);
    }
    /**
     * Add this function because we cannot simply return &(*this)[index], it is invalid to dereference a device pointer on host.
     */
    template<class Derived>
    __host__ __device__ inline LValueVector<Derived>::PtrTy LValueVector<Derived>::data_ptr(size_t index) noexcept {
        assert(index < Base::getLength() && "[Error]: Index out of range");
        return Base::getDerived().data_ptr(index);
    }

    template<class Derived>
    __host__ __device__ inline LValueVector<Derived>::ConstPtrTy LValueVector<Derived>::data_ptr(size_t index) const noexcept {
        return const_cast<This&>(*this).data_ptr(index);
    }
}
