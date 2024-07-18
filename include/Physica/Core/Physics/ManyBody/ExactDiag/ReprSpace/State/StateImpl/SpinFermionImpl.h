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
    template<unsigned int Dim, unsigned int NumSite>
    inline SpinFermion<Dim, NumSite>::SpinFermion(SpinlessType spinUp_, SpinlessType spinDown_)
            : spinUp(spinUp_), spinDown(spinDown_) {}

    template<unsigned int Dim, unsigned int NumSite>
    inline bool SpinFermion<Dim, NumSite>::operator>(const This& other) const noexcept {
        if (spinUp > other.spinUp)
            return true;
        if (spinUp == other.spinUp)
            return spinDown > other.spinDown;
        return false;
    }

    template<unsigned int Dim, unsigned int NumSite>
    inline bool SpinFermion<Dim, NumSite>::operator<(const This& other) const noexcept {
        if (spinUp < other.spinUp)
            return true;
        if (spinUp == other.spinUp)
            return spinDown < other.spinDown;
        return false;
    }

    template<unsigned int Dim, unsigned int NumSite>
    SpinFermion<Dim, NumSite> SpinFermion<Dim, NumSite>::hopUp(unsigned char from, unsigned char to) const {
        auto newSpinUp = spinUp.hop(from, to);
        const bool hopFailed = newSpinUp.isVacuum() && !spinUp.isVacuum();
        if (hopFailed)
            return SpinFermion();
        return SpinFermion(std::move(newSpinUp), spinDown);
    }

    template<unsigned int Dim, unsigned int NumSite>
    SpinFermion<Dim, NumSite> SpinFermion<Dim, NumSite>::hopDown(unsigned char from, unsigned char to) const {
        auto newSpinDown = spinDown.hop(from, to);
        const bool hopFailed = newSpinDown.isVacuum() && !spinDown.isVacuum();
        if (hopFailed)
            return SpinFermion();
        return SpinFermion(spinUp, std::move(newSpinDown));
    }

    template<unsigned int Dim, unsigned int NumSite>
    inline int SpinFermion<Dim, NumSite>::hopUpSign(unsigned char from, unsigned char to) const {
        return spinUp.hopSign(from, to);
    }

    template<unsigned int Dim, unsigned int NumSite>
    inline int SpinFermion<Dim, NumSite>::hopDownSign(unsigned char from, unsigned char to) const {
        return spinDown.hopSign(from, to);
    }

    template<unsigned int Dim, unsigned int NumSite>
    SpinFermion<Dim, NumSite> SpinFermion<Dim, NumSite>::transReduce() const {
        if constexpr (Dim != 1)
            throw NotImplementedException();
        This result = *this, temp = *this;
        for (unsigned int i = 0; i < NumSite; ++i) {
            temp <<= 1;
            if (temp.getSpinUp().getOccupyBits() < result.getSpinUp().getOccupyBits())
                result = temp;
        }
        result.spinDown = result.spinDown.transReduce(result.spinUp.calcPeriod());
        return result;
    }

    template<unsigned int Dim, unsigned int NumSite>
    inline int SpinFermion<Dim, NumSite>::lShiftSign() const {
        return spinUp.lShiftSign() * spinDown.lShiftSign();
    }

    template<unsigned int Dim, unsigned int NumSite>
    inline int SpinFermion<Dim, NumSite>::calcPeriod() const noexcept {
        if constexpr (Dim != 1)
            throw NotImplementedException();
        const int result = lcm<int, false>(spinUp.calcPeriod(), spinDown.calcPeriod());
        assert(0 < result && result <= NumSite && "[Error]: Unexpected period, this is a bug");
        return result;
    }

    template<unsigned int Dim, unsigned int NumSite>
    inline void SpinFermion<Dim, NumSite>::swap(SpinFermion& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        spinUp.swap(obj.spinUp);
        spinDown.swap(obj.spinDown);
    }

    template<unsigned int Dim, unsigned int NumSite>
    inline unsigned int SpinFermion<Dim, NumSite>::getNumDoubleOccupy() const noexcept {
        return countOnes(spinUp.getOccupyBits() & spinDown.getOccupyBits());
    }
}
