/*
 * Copyright 2021-2022 WeiBo He.
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

#include "RValueVector.h"
#include "VectorBlock.h"

namespace Physica::Core {
    /**
     * \class LValueVector is base class of vectors that can be assigned to \class LValueVector
     * and other vectors can be assigned to this class.
     * In other words, you can take the address of elements in the vector.
     */
    template<class Derived>
    class LValueVector : public RValueVector<Derived> {
    public:
        using Base = RValueVector<Derived>;
        using typename Base::ScalarType;
    public:
        ~LValueVector() = default;
        /* Operators */
        LValueVector& operator=(const LValueVector& v);
        LValueVector& operator=(LValueVector&& v) noexcept;
        template<class OtherVector>
        Derived& operator=(const RValueVector<OtherVector>& v);
        template<class AnyScalar>
        Derived& operator=(const ScalarBase<AnyScalar>& s);
        [[nodiscard]] ScalarType& operator[](size_t index) { return Base::getDerived()[index]; }
        [[nodiscard]] const ScalarType& operator[](size_t index) const { return Base::getDerived()[index]; }
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t index) const { return (*this)[index]; }
        template<class PacketType>
        void writePacket(size_t index, const PacketType packet);
        template<class PacketType>
        void writePacketPartial(size_t index, const PacketType packet);
        VectorBlock<Derived> head(size_t to) { return VectorBlock<Derived>(Base::getDerived(), 0, to); }
        const VectorBlock<Derived> head(size_t to) const { return VectorBlock<Derived>(Base::getConstCastDerived(), 0, to); }
        VectorBlock<Derived> tail(size_t from) { return VectorBlock<Derived>(Base::getDerived(), from); }
        const VectorBlock<Derived> tail(size_t from) const { return VectorBlock<Derived>(Base::getConstCastDerived(), from); }
        VectorBlock<Derived> segment(size_t from, size_t to) { return VectorBlock<Derived>(Base::getDerived(), from, to); }
        const VectorBlock<Derived> segment(size_t from, size_t to) const { return VectorBlock<Derived>(Base::getConstCastDerived(), from, to); }
        /* Getters */
        [[nodiscard]] bool isZero() const;
    protected:
        LValueVector() = default;
        LValueVector(const LValueVector&) = default;
        LValueVector(LValueVector&&) noexcept = default;
    };
}

#include "LValueVectorImpl.h"
