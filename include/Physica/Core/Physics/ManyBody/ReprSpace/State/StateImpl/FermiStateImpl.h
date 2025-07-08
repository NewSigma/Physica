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

#include "Physica/Core/Math/NumberTheory/NumberTheory.h"
#include "../FermiState.h"

namespace Physica {
    template<int Dim, int NumSite>
    FermiState<Dim, NumSite>::FermiState(Spin spinUp_, Spin spinDown_)
            : spinUp(spinUp_), spinDown(spinDown_) {}

    template<int Dim, int NumSite>
    bool FermiState<Dim, NumSite>::operator>(const This& other) const noexcept {
        if (spinUp > other.spinUp)
            return true;
        if (spinUp == other.spinUp)
            return spinDown > other.spinDown;
        return false;
    }

    template<int Dim, int NumSite>
    bool FermiState<Dim, NumSite>::operator<(const This& other) const noexcept {
        if (spinUp < other.spinUp)
            return true;
        if (spinUp == other.spinUp)
            return spinDown < other.spinDown;
        return false;
    }

    template<int Dim, int NumSite>
    auto FermiState<Dim, NumSite>::hopUp(uint8_t from, uint8_t to) const noexcept -> This {
        auto newSpinUp = spinUp.hop(from, to);
        const bool hopFailed = newSpinUp.isVacuum() && !spinUp.isVacuum();
        if (hopFailed)
            return FermiState();
        return FermiState(std::move(newSpinUp), spinDown);
    }

    template<int Dim, int NumSite>
    auto FermiState<Dim, NumSite>::hopDown(uint8_t from, uint8_t to) const noexcept -> This {
        auto newSpinDown = spinDown.hop(from, to);
        const bool hopFailed = newSpinDown.isVacuum() && !spinDown.isVacuum();
        if (hopFailed)
            return FermiState();
        return FermiState(spinUp, std::move(newSpinDown));
    }

    template<int Dim, int NumSite>
    int FermiState<Dim, NumSite>::hopUpSign(uint8_t from, uint8_t to) const noexcept {
        return spinUp.hopSign(from, to);
    }

    template<int Dim, int NumSite>
    int FermiState<Dim, NumSite>::hopDownSign(uint8_t from, uint8_t to) const noexcept {
        return spinDown.hopSign(from, to);
    }

    template<int Dim, int NumSite>
    auto FermiState<Dim, NumSite>::transReduce() const -> This {
        if constexpr (Dim != 1)
            noImpl(__func__);
        This result = *this, temp = *this;
        for (int i = 0; i < NumSite; ++i) {
            temp <<= 1;
            if (temp.getSpinUp().getOccupyBits() < result.getSpinUp().getOccupyBits())
                result = temp;
        }
        result.spinDown = result.spinDown.transReduce(result.spinUp.calcPeriod());
        return result;
    }

    template<int Dim, int NumSite>
    int FermiState<Dim, NumSite>::lShiftSign() const {
        return spinUp.lShiftSign() * spinDown.lShiftSign();
    }

    template<int Dim, int NumSite>
    int FermiState<Dim, NumSite>::calcPeriod() const noexcept {
        if constexpr (Dim != 1)
            noImpl(__func__);
        const int result = lcm<int, false>(spinUp.calcPeriod(), spinDown.calcPeriod());
        assert(0 < result && result <= NumSite && "[Error]: Unexpected period, this is a bug");
        return result;
    }

    template<int Dim, int NumSite>
    void FermiState<Dim, NumSite>::swap(FermiState& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        spinUp.swap(obj.spinUp);
        spinDown.swap(obj.spinDown);
    }

    template<int Dim, int NumSite>
    int FermiState<Dim, NumSite>::getNumDoubleOccupy() const noexcept {
        return std::popcount(spinUp.getOccupyBits() & spinDown.getOccupyBits());
    }

    template<int Dim, int NumSite>
    template<RNG R>
    FermiState<Dim, NumSite> FermiState<Dim, NumSite>::random_state() {
        return This(Spin::template random_state<R>(), Spin::template random_state<R>());
    }

    template<int Dim, int NumSite>
    template<RNG R>
    FermiState<Dim, NumSite> FermiState<Dim, NumSite>::random_state(size_t numSpinUp, size_t numSpinDown) {
        return This(Spin::template random_state<R>(numSpinUp), Spin::template random_state<R>(numSpinDown));
    }
}
