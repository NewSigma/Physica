/*
 * Copyright 2021-2026 Weibo He.
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
#include "LValueVectorImpl/LVectorView.h"

namespace Physica {
    template<class V> class RealVectorL;
    /**
     * \class LValueVector is a base class for vectors whose elements you can take the address of.
     * Any vector classes can be assigned to an LValueVector.
     */
    template<class Derived>
    class LValueVector : public RValueVector<Derived> {
        using Base = RValueVector<Derived>;
        using This = LValueVector<Derived>;
    public:
        using Base::isForwardDiff;
        using Base::isReverseDiff;
        using Base::isComplex;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using typename Base::Tr;
        using typename Base::Trv;
    public:
        ~LValueVector() = default;
        /* Operators */
        This& operator=(const This& v) = delete;
        This& operator=(This&& v) = delete;

        auto operator=(Scalar auto x) noexcept -> Derived&;
        void operator+=(Scalar auto x) noexcept;
        void operator-=(Scalar auto x) noexcept;
        void operator*=(Scalar auto x) noexcept;
        void operator/=(Scalar auto x) noexcept;

        template<ExecutePolicy P = Sequential>
        auto operator=(const Vector auto& v) -> Derived&;
        void operator+=(const Vector auto& v);
        void operator-=(const Vector auto& v);

        [[nodiscard]] decltype(auto) operator[](this auto&&, size_t index);
        /* Operations */
        [[nodiscard]] decltype(auto) calc(size_t index) const;
        [[nodiscard]] Tv calc_value(size_t index) const;
        void writePacket(Packet auto packet, size_t index) noexcept;
        void writePacket(Packet auto packet, size_t index, size_t count) noexcept;
        [[nodiscard]] constexpr auto view(this auto&&) noexcept;

        [[nodiscard]] CoDiff<T> sum() const;

        void toNextMean(size_t lastNumSample, const Vector auto& sample) noexcept;
        void toNextVariance(Derived& mean, size_t lastNumSample, const Vector auto& sample) noexcept;
        void reverse(const auto& grad) const noexcept;

        template<size_t Length = Dynamic>
        [[nodiscard]] auto head(this auto&&, size_t to = Length) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] auto tail(this auto&&, size_t from) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] auto segment(this auto&&, size_t from, size_t to) noexcept;

        [[nodiscard]] auto reshape_like(this auto&& self, const Matrix auto& mat) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto reshape_col(this auto&& self, size_t row, size_t col) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto reshape_row(this auto&& self, size_t row, size_t col) noexcept;

        void zeros() noexcept;
        void clamp_min(Tv minimum);
        void clamp_max(Tv maximum);
        void toUnit();
        void standardize();
        using Base::householder;
        Tr householder();
        template<RNG R>
        void random_uniform();
        template<RNG R>
        void random_normal();
        template<RNG R>
        void random_any(auto& distribution);

        [[nodiscard]] decltype(auto) reals() noexcept;
        [[nodiscard]] decltype(auto) reals() const noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] auto grads() const noexcept;
        /* Getters */
        [[nodiscard]] auto data_ptr(this auto&&, size_t index) noexcept;
        [[nodiscard]] auto& front(this auto&&) noexcept;
        [[nodiscard]] auto& back(this auto&&) noexcept;
    protected:
        LValueVector() = default;
        LValueVector(const This&) = default;
        LValueVector(This&&) noexcept = default;
    };
}

#include "LValueVectorImpl/LValueVectorImpl.h"
#include "LValueVectorImpl/VectorConvert/RealVector.h"
#include "LValueVectorImpl/Sincos.h"
