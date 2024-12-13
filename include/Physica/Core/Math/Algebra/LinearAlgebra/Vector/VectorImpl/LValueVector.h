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
        using Base::isReverseDiff;
    protected:
        using PtrTy = ScalarType::PtrTy;
        using ConstPtrTy = ScalarType::ConstPtrTy;
        using RefTy = ScalarType::RefTy;
        using ConstRefTy = ScalarType::ConstRefTy;
    public:
        ~LValueVector() = default;
        /* Operators */
        inline This& operator=(const This& v);
        inline This& operator=(This&& v);

        template<Scalar T> inline Derived& operator=(const T& x);
        template<Scalar T> void operator+=(const T& x) { Base::getDerived() = Base::getDerived() + x; }
        template<Scalar T> void operator-=(const T& x) { Base::getDerived() = Base::getDerived() - x; }
        template<Scalar T> void operator*=(const T& x) { Base::getDerived() = Base::getDerived() * x; }
        template<Scalar T> void operator/=(const T& x) { Base::getDerived() = Base::getDerived() / x; }

        template<Vector T, class Executor = SequentialExecutor>
        inline Derived& operator=(const T& v_);
        template<Vector T> void operator+=(const T& v) { Base::getDerived() = Base::getDerived() + v; }
        template<Vector T> void operator-=(const T& v) { Base::getDerived() += -v; }
        [[nodiscard]] inline RefTy operator[](size_t index);
        [[nodiscard]] inline ConstRefTy operator[](size_t index) const;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t index) const;
        template<class AnyPacket> void writePacket(size_t index, const AnyPacket packet);
        template<class AnyPacket> void writePacketPartial(size_t index, size_t count, const AnyPacket packet);

        [[nodiscard]] CoDiff<ScalarType> sum() const;

        template<Vector T>
        void reverse(const T& v) const noexcept requires(isReverseDiff);

        template<size_t Length = Dynamic>
        [[nodiscard]] inline auto head(size_t to) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline const auto head(size_t to) const noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline auto tail(size_t from) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline const auto tail(size_t from) const noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline auto segment(size_t from, size_t to) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline const auto segment(size_t from, size_t to) const noexcept;

        inline void toUnit();
        template<RandomGenerator R>
        inline void random_uniform();
        template<RandomGenerator R>
        inline void random_normal();
        template<class Distribution, RandomGenerator R>
        inline void random_any(Distribution& dist);
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
