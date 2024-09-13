/*
 * Copyright 2021-2024 Weibo He.
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
#include "LValueVectorImpl/LVectorBlock.h"

namespace Physica::Core {
    /**
     * \class LValueVector is a base class for vectors. You can take the address of elements in an LValueVector.
     * An LValueVector can be assigned to an LValueVector, and other vector classes can be assigned to an LValueVector.
     */
    template<class Derived>
    class LValueVector : public RValueVector<Derived> {
        using Base = RValueVector<Derived>;
        template<size_t Length>
        using BlockType = LVectorBlock<Derived, Length>;
    public:
        
        using typename Base::ScalarType;
    public:
        ~LValueVector() = default;
        /* Operators */
        inline LValueVector& operator=(const LValueVector& v);
        inline LValueVector& operator=(LValueVector&& v) noexcept;
        template<class OtherVector, class Executor = SequentialExecutor>
        inline Derived& operator=(const RValueVector<OtherVector>& v_);
        inline Derived& operator=(const ScalarType& s);
        [[nodiscard]] ScalarType& operator[](size_t index) { return *data_ptr(index); }
        [[nodiscard]] const ScalarType& operator[](size_t index) const { return *data_ptr(index); }
        void operator+=(const ScalarType& s) { (*this) = (*this) + s; }
        void operator-=(const ScalarType& s) { (*this) = (*this) - s; }
        void operator*=(const ScalarType& s) { (*this) = (*this) * s; }
        void operator/=(const ScalarType& s) { (*this) = (*this) / s; }
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t index) const { return *data_ptr(index); }
        template<class AnyPacket> void writePacket(size_t index, const AnyPacket packet);
        template<class AnyPacket> void writePacketPartial(size_t index, size_t count, const AnyPacket packet);

        template<size_t Length = Dynamic>
        [[nodiscard]] inline BlockType<Length> head(size_t to);
        template<size_t Length = Dynamic>
        [[nodiscard]] inline const BlockType<Length> head(size_t to) const;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline BlockType<Length> tail(size_t from);
        template<size_t Length = Dynamic>
        [[nodiscard]] inline const BlockType<Length> tail(size_t from) const;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline BlockType<Length> segment(size_t from, size_t to);
        template<size_t Length = Dynamic>
        [[nodiscard]] inline const BlockType<Length> segment(size_t from, size_t to) const;

        inline void toUnit();
        template<class RandomGenerator>
        inline void random_uniform(RandomGenerator& gen);
        template<class RandomGenerator>
        inline void random_normal(RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        inline void random_any(Distribution& dist, RandomGenerator& gen);
        /* Getters */
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

#include "LValueVectorImpl/LValueVectorImpl.h"
#include "Sincos.h"
