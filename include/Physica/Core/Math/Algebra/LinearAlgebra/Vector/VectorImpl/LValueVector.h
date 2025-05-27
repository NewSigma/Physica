/*
 * Copyright 2021-2025 Weibo He.
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

namespace Physica {
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
        using typename Base::Tv;
        using typename Base::Tr;
        using PtrTy = ScalarType::PtrTy;
        using ConstPtrTy = ScalarType::ConstPtrTy;
        using RefTy = ScalarType::RefTy;
        using ConstRefTy = ScalarType::ConstRefTy;
    public:
        ~LValueVector() = default;
        /* Operators */
        inline This& operator=(const This& v);
        inline This& operator=(This&& v);

        template<Scalar T> inline Derived& operator=(const T& x) requires(!isReverseDiff || !ReverseDiff<T>);
        template<Scalar T> void operator+=(const T& x) { Base::getDerived() = Base::getDerived() + x; }
        template<Scalar T> void operator-=(const T& x) { Base::getDerived() = Base::getDerived() - x; }
        template<Scalar T> void operator*=(const T& x) { Base::getDerived() = Base::getDerived() * x; }
        template<Scalar T> void operator/=(const T& x) { Base::getDerived() = Base::getDerived() / x; }

        template<Vector V, ExecutePolicy P = Sequential>
        inline Derived& operator=(const V& v_) requires(!CUDA<V>);
        template<Vector V> inline void operator+=(const V& v);
        template<Vector V> inline void operator-=(const V& v);

        [[nodiscard]] inline RefTy operator[](size_t index);
        [[nodiscard]] inline ConstRefTy operator[](size_t index) const;
        /* Operations */
        [[nodiscard]] ConstRefTy calc(size_t index) const;
        [[nodiscard]] Tv calc_value(size_t index) const;
        template<Packet Pack> void writePacket(size_t index, const Pack packet);
        template<Packet Pack> void writePacketPartial(size_t index, size_t count, const Pack packet);

        [[nodiscard]] CoDiff<ScalarType> sum() const;

        template<class T>
        void reverse(const T& grad) const noexcept requires(isReverseDiff);

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

        template<Matrix M>
        [[nodiscard]] auto reshape(const M& mat) noexcept;
        template<Matrix M>
        [[nodiscard]] const auto reshape(const M& mat) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto reshape_col(size_t row, size_t col) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] const auto reshape_col(size_t row, size_t col) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto reshape_row(size_t row, size_t col) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] const auto reshape_row(size_t row, size_t col) const noexcept;

        void clamp_min(const Tv& minimum);
        void clamp_max(const Tv& maximum);
        void toUnit();
        using Base::householder;
        Tr householder();

        template<RNG R>
        inline void random_uniform();
        template<RNG R>
        inline void random_normal();
        template<RNG R, class Distribution>
        inline void random_any(Distribution& dist);

        template<int GradOrder = 1>
        auto grads() const noexcept;
        /* Getters */
        [[nodiscard]] inline PtrTy data_ptr(size_t index) noexcept;
        [[nodiscard]] inline ConstPtrTy data_ptr(size_t index) const noexcept;
    protected:
        LValueVector() = default;
        LValueVector(const This&) = default;
        LValueVector(This&&) noexcept = default;
    };
}

#include "LValueVectorImpl/LValueVectorImpl.h"
#include "LValueVectorImpl/Sincos.h"
