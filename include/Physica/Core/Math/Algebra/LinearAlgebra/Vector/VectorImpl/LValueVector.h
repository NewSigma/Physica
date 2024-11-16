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
        using This = LValueVector<Derived>;
        template<size_t Length>
        using BlockType = LVectorBlock<Derived, Length>;
    public:
        using typename Base::ScalarType;
        using Base::isForwardDiff;
    protected:
        using PtrTy = typename ScalarType::PtrTy;
        using ConstPtrTy = typename ScalarType::ConstPtrTy;
        using RefTy = typename ScalarType::RefTy;
        using ConstRefTy = typename ScalarType::ConstRefTy;
    public:
        ~LValueVector() = default;
        /* Operators */
        inline This& operator=(const This& v);
        inline This& operator=(This&& v);

        template<Scalar T>
        inline Derived& operator=(const T& x);
        void operator+=(const ScalarType& s) { (*this) = (*this) + s; }
        void operator-=(const ScalarType& s) { (*this) = (*this) - s; }
        void operator*=(const ScalarType& s) { (*this) = (*this) * s; }
        void operator/=(const ScalarType& s) { (*this) = (*this) / s; }

        template<class OtherVector, class Executor = SequentialExecutor>
        inline Derived& operator=(const RValueVector<OtherVector>& v_);
        template<class OtherVector> void operator+=(const RValueVector<OtherVector>& v) { (*this) = (*this) + v; }
        template<class OtherVector> void operator-=(const RValueVector<OtherVector>& v) { Base::getDerived() += (-v.getDerived()); }
        [[nodiscard]] inline RefTy operator[](size_t index);
        [[nodiscard]] inline ConstRefTy operator[](size_t index) const;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t index) const;
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
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(size_t index) noexcept;
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(size_t index) const noexcept;
    protected:
        LValueVector() = default;
        LValueVector(const This&) = default;
        LValueVector(This&&) noexcept = default;
    };
}

#include "LValueVectorImpl/LValueVectorImpl.h"
#include "Sincos.h"
