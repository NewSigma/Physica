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

namespace Physica::Core {
    template<class Derived>
    inline LValueVector<Derived>& LValueVector<Derived>::operator=(const LValueVector& v) {
        return Base::getDerived() = v.getDerived();
    }

    template<class Derived>
    inline LValueVector<Derived>& LValueVector<Derived>::operator=(LValueVector&& v) noexcept {
        return Base::getDerived() = std::move(v.getDerived());
    }

    template<class Derived>
    template<class OtherVector, class Executor>
    inline Derived& LValueVector<Derived>::operator=(const RValueVector<OtherVector>& v_) {
        Derived& v1 = Base::getDerived();
        const auto& v = v_.getDerived();
        v1.resize(v.getLength());
        v.template assignTo<Derived, Executor>(v1);
        return v1;
    }

    template<class Derived>
    template<class OtherScalar>
    inline Derived& LValueVector<Derived>::operator=(const ScalarBase<OtherScalar>& s) {
        const auto x = ScalarType(s.getDerived());
        for (size_t i = 0; i < Base::getLength(); ++i)
            (*this)[i] = x;
        return Base::getDerived();
    }

    template<class Derived>
    inline typename LValueVector<Derived>::RefTy LValueVector<Derived>::operator[](size_t index) {
        return *data_ptr(index);
    }

    template<class Derived>
    inline typename LValueVector<Derived>::ConstRefTy LValueVector<Derived>::operator[](size_t index) const {
        return const_cast<This&>(*this).operator[](index);
    }

    template<class Derived>
    typename LValueVector<Derived>::ScalarType LValueVector<Derived>::calc(size_t index) const {
        return ScalarType(*data_ptr(index));
    }

    template<class Derived>
    template<class AnyPacket>
    void LValueVector<Derived>::writePacket(size_t index, const AnyPacket packet) {
        ScalarType buffer[AnyPacket::size()];
        packet.store(buffer);
        for (size_t i = 0; i < AnyPacket::size(); ++i, ++index)
            (*this)[index] = buffer[i];
    }

    template<class Derived>
    template<class AnyPacket>
    void LValueVector<Derived>::writePacketPartial(size_t index, size_t count, const AnyPacket packet) {
        ScalarType buffer[AnyPacket::size()];
        packet.store(buffer);
        for (size_t i = 0; i < count; ++i, ++index)
            (*this)[index] = buffer[i];
    }

    template<class Derived>
    template<size_t Length>
    inline LVectorBlock<Derived, Length> LValueVector<Derived>::head(size_t to) {
        return {Base::getDerived(), 0, to};
    }

    template<class Derived>
    template<size_t Length>
    inline const LVectorBlock<Derived, Length> LValueVector<Derived>::head(size_t to) const {
        return {Base::getConstCastDerived(), 0, to};
    }

    template<class Derived>
    template<size_t Length>
    inline LVectorBlock<Derived, Length> LValueVector<Derived>::tail(size_t from) {
        return {Base::getDerived(), from};
    }

    template<class Derived>
    template<size_t Length>
    inline const LVectorBlock<Derived, Length> LValueVector<Derived>::tail(size_t from) const {
        return {Base::getConstCastDerived(), from};
    }

    template<class Derived>
    template<size_t Length>
    inline LVectorBlock<Derived, Length> LValueVector<Derived>::segment(size_t from, size_t to) {
        return {Base::getDerived(), from, to};
    }

    template<class Derived>
    template<size_t Length>
    inline const LVectorBlock<Derived, Length> LValueVector<Derived>::segment(size_t from, size_t to) const {
        return {Base::getConstCastDerived(), from, to};
    }

    template<class Derived>
    inline void LValueVector<Derived>::toUnit() {
        Base::getDerived() *= reciprocal(Base::getDerived().norm());
    }

    template<class Derived>
    template<class RandomGenerator>
    inline void LValueVector<Derived>::random_uniform(RandomGenerator& gen) {
        for (size_t i = 0; i < this->getLength(); ++i)
            this->operator[](i) = ScalarType::random_uniform(gen);
    }

    template<class Derived>
    template<class RandomGenerator>
    inline void LValueVector<Derived>::random_normal(RandomGenerator& gen) {
        for (size_t i = 0; i < this->getLength(); ++i)
            this->operator[](i) = ScalarType::random_normal(gen);
    }

    template<class Derived>
    template<class Distribution, class RandomGenerator>
    inline void LValueVector<Derived>::random_any(Distribution& dist, RandomGenerator& gen) {
        for (size_t i = 0; i < this->getLength(); ++i)
            this->operator[](i) = ScalarType::random_any(dist, gen);
    }
    /**
     * Add this function because we cannot simply return &(*this)[index], it is invalid to dereference a device pointer on host.
     */
    template<class Derived>
    __host__ __device__ inline typename LValueVector<Derived>::PtrTy LValueVector<Derived>::data_ptr(size_t index) noexcept {
        return Base::getDerived().data_ptr(index);
    }

    template<class Derived>
    __host__ __device__ inline typename LValueVector<Derived>::ConstPtrTy LValueVector<Derived>::data_ptr(size_t index) const noexcept {
        return const_cast<This&>(*this).data_ptr(index);
    }

    template<class Derived, class OtherDerived>
    inline void operator+=(LValueVector<Derived>& v1, const RValueVector<OtherDerived>& v2) {
        v1 = v1 + v2;
    }

    template<class Derived, class OtherDerived>
    inline void operator-=(LValueVector<Derived>& v1, const RValueVector<OtherDerived>& v2) {
        v1.getDerived() += (-v2.getDerived());
    }
}
