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
        using Base::isForwardDiff;
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using typename Base::Tr;
        using typename Base::Trv;
        using PtrTy = T::PtrTy;
        using ConstPtrTy = T::ConstPtrTy;
        using RefTy = T::RefTy;
        using ConstRefTy = T::ConstRefTy;
    public:
        ~LValueVector() = default;
        /* Operators */
        This& operator=(const This& v) = delete;
        This& operator=(This&& v) = delete;

        Derived& operator=(const Scalar auto& x);
        void operator+=(const Scalar auto& x) { Base::getDerived() = Base::getDerived() + x; }
        void operator-=(const Scalar auto& x) { Base::getDerived() = Base::getDerived() - x; }
        void operator*=(const Scalar auto& x) { Base::getDerived() = Base::getDerived() * x; }
        void operator/=(const Scalar auto& x) { Base::getDerived() = Base::getDerived() / x; }

        template<ExecutePolicy P = Sequential>
        Derived& operator=(const Vector auto& v);
        void operator+=(const Vector auto& v);
        void operator-=(const Vector auto& v);

        [[nodiscard]] RefTy operator[](size_t index);
        [[nodiscard]] ConstRefTy operator[](size_t index) const;
        /* Operations */
        [[nodiscard]] ConstRefTy calc(size_t index) const;
        [[nodiscard]] Tv calc_value(size_t index) const;
        void writePacket(size_t index, Packet auto packet);
        void writePacketPartial(size_t index, size_t count, Packet auto packet);

        [[nodiscard]] CoDiff<T> sum() const;

        void toNextMean(size_t lastNumSample, const Vector auto& sample) noexcept;
        void toNextVariance(Derived& mean, size_t lastNumSample, const Vector auto& sample) noexcept;
        void reverse(const auto& grad) const noexcept;

        template<size_t Length = Dynamic>
        [[nodiscard]] auto head(size_t to) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] const auto head(size_t to) const noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] auto tail(size_t from) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] const auto tail(size_t from) const noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] auto segment(size_t from, size_t to) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] const auto segment(size_t from, size_t to) const noexcept;

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
        void standardize();
        using Base::householder;
        Tr householder();

        void zeros() noexcept;
        template<RNG R = Random<>>
        void random_uniform();
        template<RNG R = Random<>>
        void random_normal();
        template<RNG R = Random<>>
        void random_any(auto& distribution);

        template<int GradOrder = 1>
        auto grads() const noexcept;
        /* Getters */
        [[nodiscard]] PtrTy data_ptr(size_t index) noexcept;
        [[nodiscard]] ConstPtrTy data_ptr(size_t index) const noexcept;
    protected:
        LValueVector() = default;
        LValueVector(const This&) = default;
        LValueVector(This&&) noexcept = default;
    };
}

#include "LValueVectorImpl/LValueVectorImpl.h"
#include "LValueVectorImpl/Sincos.h"
