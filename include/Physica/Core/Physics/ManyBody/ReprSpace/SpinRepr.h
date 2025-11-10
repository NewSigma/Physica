/*
 * Copyright 2024 Weibo He.
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

#include "State/SpinState.h"
#include "ReprBasis.h"

namespace Physica {
    template<int Dim, int NumSite>
    class SpinRepr : public ReprBasis<SpinRepr<Dim, NumSite>> {
        using This = SpinRepr<Dim, NumSite>;
        using Base = ReprBasis<This>;
    public:
        using typename Base::StateType;
    public:
        SpinRepr() = default;
        SpinRepr(const This&) = default;
        SpinRepr(This&&) noexcept = default;
        ~SpinRepr() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] StateType operator[](size_t index) const noexcept;
        [[nodiscard]] size_t operator[](StateType state) const noexcept;
        /* Operations */
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] constexpr static int getNumSpin() noexcept { return NumSite; }
        [[nodiscard]] constexpr static size_t getNumState() noexcept { return size_t(1) << NumSite; }
    };

    template<int Dim, int NumSite>
    auto SpinRepr<Dim, NumSite>::operator[](size_t index) const noexcept -> StateType {
        assert(index < getNumState() && "[Error]: Index out of range");
        return StateType(index);
    }

    template<int Dim, int NumSite>
    size_t SpinRepr<Dim, NumSite>::operator[](StateType state) const noexcept {
        assert(state.getOccupyBits() < getNumState() && "[Error]: Index out of range");
        return state.getOccupyBits();
    }

    template<int Dim, int NumSite>
    void SpinRepr<Dim, NumSite>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
    }
}

namespace Physica {
    template<int Dim_, int NumSite>
    class Traits<SpinRepr<Dim_, NumSite>> {
    public:
        using StateType = SpinState<Dim_, NumSite>;
        constexpr static int Dim = Dim_;
        constexpr static bool IsTransInvariant = false;
    };
}
