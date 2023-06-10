/*
 * Copyright 2022-2023 WeiBo He.
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
    template<class OtherVector>
    inline Derived& LValueVector<Derived>::operator=(const RValueVector<OtherVector>& v) {
        Base::getDerived().resize(v.getLength());
        v.getDerived().assignTo(Base::getDerived());
        return Base::getDerived();
    }

    template<class Derived>
    template<class AnyScalar>
    inline Derived& LValueVector<Derived>::operator=(const ScalarBase<AnyScalar>& s) {
        for (size_t i = 0; i < Base::getLength(); ++i)
            (*this)[i] = s.getDerived();
        return Base::getDerived();
    }

    template<class Derived>
    template<class PacketType>
    void LValueVector<Derived>::writePacket(size_t index, const PacketType packet) {
        using TrivialType = typename ScalarType::TrivialType;
        TrivialType buffer[PacketType::size()];
        packet.store(buffer);
        for (int i = 0; i < PacketType::size(); ++i, ++index)
            (*this)[index] = buffer[i];
    }

    template<class Derived>
    template<class PacketType>
    void LValueVector<Derived>::writePacketPartial(size_t index, size_t count, const PacketType packet) {
        using TrivialType = typename ScalarType::TrivialType;
        TrivialType buffer[PacketType::size()];
        packet.store(buffer);
        for (size_t i = 0; i < count; ++i, ++index)
            (*this)[index] = buffer[i];
    }

    template<class Derived>
    void LValueVector<Derived>::toUnit() {
        Base::getDerived() *= reciprocal(Base::getDerived().norm());
    }

    template<class Derived>
    template<class RandomGenerator>
    void LValueVector<Derived>::random_uniform(RandomGenerator& gen) {
        for (size_t i = 0; i < this->getLength(); ++i)
            this->operator[](i) = ScalarType::random_uniform(gen);
    }

    template<class Derived>
    template<class RandomGenerator>
    void LValueVector<Derived>::random_normal(RandomGenerator& gen) {
        for (size_t i = 0; i < this->getLength(); ++i)
            this->operator[](i) = ScalarType::random_normal(gen);
    }

    template<class Derived>
    bool LValueVector<Derived>::isZero() const {
        for(size_t i = 0; i < Base::getLength(); ++i)
            if (!(*this)[i].isZero())
                return false;
        return true;
    }

    template<class Derived, class OtherDerived>
    inline void operator+=(LValueVector<Derived>& v1, const RValueVector<OtherDerived>& v2) {
        v1 = v1 + v2;
    }

    template<class Derived, class OtherDerived>
    inline void operator-=(LValueVector<Derived>& v1, const RValueVector<OtherDerived>& v2) {
        v1.getDerived() += (-v2.getDerived());
    }

    template<class VectorType, class ScalarType>
    inline void operator+=(LValueVector<VectorType>& v, const ScalarBase<ScalarType>& n) { v = v + n; }

    template<class VectorType, class ScalarType>
    inline void operator-=(LValueVector<VectorType>& v, const ScalarBase<ScalarType>& n) { v = v - n; }

    template<class VectorType, class ScalarType>
    inline void operator*=(LValueVector<VectorType>& v, const ScalarBase<ScalarType>& n) { v = v * n; }

    template<class VectorType, class ScalarType>
    inline void operator/=(LValueVector<VectorType>& v, const ScalarBase<ScalarType>& n) { v = v / n; }
}
