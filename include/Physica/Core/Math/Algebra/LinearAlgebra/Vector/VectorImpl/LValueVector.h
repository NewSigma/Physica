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
#include "LVectorBlock.h"

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
        inline LValueVector& operator=(const LValueVector& v);
        inline LValueVector& operator=(LValueVector&& v) noexcept;
        template<class OtherVector>
        inline Derived& operator=(const RValueVector<OtherVector>& v);
        template<class AnyScalar>
        inline Derived& operator=(const ScalarBase<AnyScalar>& s);
        [[nodiscard]] ScalarType& operator[](size_t index) { return *data_ptr(index); }
        [[nodiscard]] const ScalarType& operator[](size_t index) const { return *data_ptr(index); }
        template<class T> void operator+=(const ScalarBase<T>& s) { (*this) = (*this) + s.getDerived(); }
        template<class T> void operator-=(const ScalarBase<T>& s) { (*this) = (*this) - s.getDerived(); }
        template<class T> void operator*=(const ScalarBase<T>& s) { (*this) = (*this) * s.getDerived(); }
        template<class T> void operator/=(const ScalarBase<T>& s) { (*this) = (*this) / s.getDerived(); }
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t index) const { return *data_ptr(index); }
        template<class PacketType> void writePacket(size_t index, const PacketType packet);
        template<class PacketType> void writePacketPartial(size_t index, size_t count, const PacketType packet);

        [[nodiscard]] inline LVectorBlock<Derived> head(size_t to);
        [[nodiscard]] inline const LVectorBlock<Derived> head(size_t to) const;
        [[nodiscard]] inline LVectorBlock<Derived> tail(size_t from);
        [[nodiscard]] inline const LVectorBlock<Derived> tail(size_t from) const;
        [[nodiscard]] inline LVectorBlock<Derived> segment(size_t from, size_t to);
        [[nodiscard]] inline const LVectorBlock<Derived> segment(size_t from, size_t to) const;

        inline void toUnit();
        template<class RandomGenerator>
        inline void random_uniform(RandomGenerator& gen);
        template<class RandomGenerator>
        inline void random_normal(RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        inline void random_any(Distribution& dist, RandomGenerator& gen);
        /* Getters */
        [[nodiscard]] bool isZero() const;
        [[nodiscard]] __host__ __device__ inline ScalarType* data_ptr(size_t index);
        [[nodiscard]] __host__ __device__ inline const ScalarType* data_ptr(size_t index) const;
    protected:
        LValueVector() = default;
        LValueVector(const LValueVector&) = default;
        LValueVector(LValueVector&&) noexcept = default;
    };

    template<class Derived, class OtherDerived>
    inline void operator+=(LValueVector<Derived>& v1, const RValueVector<OtherDerived>& v2);

    template<class Derived, class OtherDerived>
    inline void operator-=(LValueVector<Derived>& v1, const RValueVector<OtherDerived>& v2);
}

#include "LValueVectorImpl.h"
#include "Sincos.h"
