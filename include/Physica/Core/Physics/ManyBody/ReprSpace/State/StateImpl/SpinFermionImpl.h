/*
 * Copyright 2024 WeiBo He.
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

namespace Physica::Core {
    template<unsigned int Dim>
    inline SpinFermion<Dim>::SpinFermion(SpinlessType spinUp_, SpinlessType spinDown_)
            : spinUp(spinUp_), spinDown(spinDown_) {}

    template<unsigned int Dim>
    inline bool SpinFermion<Dim>::operator>(const This& other) const noexcept {
        if (spinUp > other.spinUp)
            return true;
        if (spinUp == other.spinUp)
            return spinDown > other.spinDown;
        return false;
    }

    template<unsigned int Dim>
    inline bool SpinFermion<Dim>::operator<(const This& other) const noexcept {
        if (spinUp < other.spinUp)
            return true;
        if (spinUp == other.spinUp)
            return spinDown < other.spinDown;
        return false;
    }

    template<unsigned int Dim>
    SpinFermion<Dim> SpinFermion<Dim>::hopUp(unsigned char from, unsigned char to) const {
        auto newSpinUp = spinUp.hop(from, to);
        const bool hopFailed = newSpinUp.isVacuum() && !spinUp.isVacuum();
        if (hopFailed)
            return SpinFermion();
        return SpinFermion(std::move(newSpinUp), spinDown);
    }

    template<unsigned int Dim>
    SpinFermion<Dim> SpinFermion<Dim>::hopDown(unsigned char from, unsigned char to) const {
        auto newSpinDown = spinDown.hop(from, to);
        const bool hopFailed = newSpinDown.isVacuum() && !spinDown.isVacuum();
        if (hopFailed)
            return SpinFermion();
        return SpinFermion(spinUp, std::move(newSpinDown));
    }

    template<unsigned int Dim>
    SpinFermion<Dim> SpinFermion<Dim>::transReduce() const {
        const int numSite = getNumSite();
        This result = *this, temp = *this;
        for (int i = 0; i < numSite; ++i) {
            temp <<= 1;
            if (temp.getSpinUp().getOccupyBits() < result.getSpinUp().getOccupyBits())
                result = temp;
        }
        result.spinDown = result.spinDown.transReduce(result.spinUp.calcPeriod());
        return result;
    }

    template<unsigned int Dim>
    inline int SpinFermion<Dim>::calcPeriod() const noexcept {
        const int result = lcm<int, false>(spinUp.calcPeriod(), spinDown.calcPeriod());
        assert(0 < result && result <= getNumSite() && "[Error]: Unexpected period, this is a bug");
        return result;
    }

    template<unsigned int Dim>
    inline void SpinFermion<Dim>::swap(SpinFermion& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        spinUp.swap(obj.spinUp);
        spinDown.swap(obj.spinDown);
    }

    template<unsigned int Dim>
    inline unsigned int SpinFermion<Dim>::getNumPairedElectron() const noexcept {
        return countOnes(spinUp.getOccupyBits() & spinDown.getOccupyBits());
    }
}
