/*
 * Copyright 2022 WeiBo He.
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

#include "LValueVector.h"

namespace Physica::Core {
    /**
     * \class ContinuousVector has its element on memory continuously
     */
    template<class Derived>
    class ContinuousVector : public LValueVector<Derived> {
    public:
        using Base = LValueVector<Derived>;
        using typename Base::ScalarType;
    private:
        using TrivialType = typename ScalarType::TrivialType;
        using ConstPointerType = typename std::add_pointer<typename std::add_const<TrivialType>::type>::type;
        using PointerType = typename std::add_pointer<TrivialType>::type;
    public:
        ~ContinuousVector() = default;
        /* Operators */
        ContinuousVector& operator=(const ContinuousVector& v);
        ContinuousVector& operator=(ContinuousVector&& v) noexcept;
        using Base::operator=;
        /* Operations */
        template<class PacketType>
        [[nodiscard]] PacketType packet(size_t index) const;
        template<class PacketType>
        [[nodiscard]] PacketType packetPartial(size_t index) const;
        template<class PacketType>
        void writePacket(size_t index, const PacketType packet);
        template<class PacketType>
        void writePacketPartial(size_t index, const PacketType packet);
    protected:
        ContinuousVector() = default;
        ContinuousVector(const ContinuousVector&) = default;
        ContinuousVector(ContinuousVector&&) noexcept = default;
    private:
        PointerType data_ptr(size_t index) { return reinterpret_cast<PointerType>(&(*this)[index]); }
        ConstPointerType data_ptr(size_t index) const { return reinterpret_cast<ConstPointerType>(&(*this)[index]); }
    };

    template<class Derived>
    ContinuousVector<Derived>& ContinuousVector<Derived>::operator=(const ContinuousVector<Derived>& v) {
        Base::operator=(v);
        return *this;
    }

    template<class Derived>
    ContinuousVector<Derived>& ContinuousVector<Derived>::operator=(ContinuousVector<Derived>&& v) noexcept {
        Base::operator=(std::forward<Base>(v));
        return *this;
    }

    template<class Derived>
    template<class PacketType>
    PacketType ContinuousVector<Derived>::packet(size_t index) const {
        assert(index + PacketType::size() <= Base::getLength());
        if constexpr (std::is_same_v<TrivialType, typename Internal::Traits<PacketType>::ScalarType>) {
            PacketType packet{};
            packet.load(data_ptr(index));
            return packet;
        }
        else
            return Base::template packet<PacketType>(index);
    }

    template<class Derived>
    template<class PacketType>
    PacketType ContinuousVector<Derived>::packetPartial(size_t index) const {
        if constexpr (std::is_same_v<TrivialType, typename Internal::Traits<PacketType>::ScalarType>) {
            PacketType packet{};
            const size_t count = Base::getLength() - index;
            packet.load_partial(count, data_ptr(index));
            return packet;
        }
        else
            return Base::template packetPartial<PacketType>(index);
    }

    template<class Derived>
    template<class PacketType>
    void ContinuousVector<Derived>::writePacket(size_t index, const PacketType packet) {
        if constexpr (std::is_same_v<TrivialType, typename Internal::Traits<PacketType>::ScalarType>)
            packet.store(data_ptr(index));
        else
            Base::template writePacket(index, packet);
    }

    template<class Derived>
    template<class PacketType>
    void ContinuousVector<Derived>::writePacketPartial(size_t index, const PacketType packet) {
        if constexpr (std::is_same_v<TrivialType, typename Internal::Traits<PacketType>::ScalarType>)
            packet.store_partial(Base::getLength() - index, data_ptr(index));
        else
            Base::template writePacketPartial(index, packet);
    }
}
