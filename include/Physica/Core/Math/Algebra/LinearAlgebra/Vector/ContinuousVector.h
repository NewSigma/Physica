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
        /* Friends */
        template<class T>
        friend std::ostream& operator<<(std::ostream& os, const ContinuousVector<T>& v);
        template<class T>
        friend std::istream& operator>>(std::istream& is, ContinuousVector<T>& v);
    };

    template<class Derived>
    std::ostream& operator<<(std::ostream& os, const ContinuousVector<Derived>& v);

    template<class Derived>
    std::istream& operator>>(std::istream& is, ContinuousVector<Derived>& v);

    template<class Derived, class OtherDerived>
    void operator+=(ContinuousVector<Derived>& v1, const RValueVector<OtherDerived>& v2);
}

#include "ContinuousVectorImpl.h"
