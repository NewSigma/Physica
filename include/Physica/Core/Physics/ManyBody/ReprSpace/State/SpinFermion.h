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

#include "SpinlessFermion.h"
#include "Physica/Core/Math/NumberTheory/NumberTheory.h"

namespace Physica::Core {
    class SpinFermion;

    namespace Internal {
        template<>
        class Traits<SpinFermion> {
        public:
            constexpr static unsigned int SiteDOF = 4;
        };
    }

    class PHYSICA_API SpinFermion : public State<SpinFermion> {
        using This = SpinFermion;

        SpinlessFermion spinUp;
        SpinlessFermion spinDown;
    public:
        SpinFermion() = default;
        inline SpinFermion(SpinlessFermion spinUp_, SpinlessFermion spinDown_);
        SpinFermion(const This&) = default;
        SpinFermion(This&&) noexcept = default;
        ~SpinFermion() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] bool operator==(const This& other) const noexcept { return spinUp == other.spinUp && spinDown == other.spinDown; }
        [[nodiscard]] bool operator!=(const This& other) const noexcept { return !(*this == other); }
        [[nodiscard]] inline bool operator>(const This& other) const noexcept;
        [[nodiscard]] inline bool operator<(const This& other) const noexcept;
        [[nodiscard]] bool operator>=(const This& other) const noexcept { return !((*this) < other); }
        [[nodiscard]] bool operator<=(const This& other) const noexcept { return !((*this) > other); }
        [[nodiscard]] This operator<<(int shift) const noexcept { return SpinFermion(spinUp << shift, spinDown << shift); }
        [[nodiscard]] This operator>>(int shift) const noexcept { return SpinFermion(spinUp >> shift, spinDown >> shift); }
        void operator<<=(int shift) noexcept { (*this) = (*this) << shift; }
        void operator>>=(int shift) noexcept { (*this) = (*this) >> shift; }
        /* Operations */
        [[nodiscard]] SpinFermion hopUp(unsigned char from, unsigned char to) const;
        [[nodiscard]] SpinFermion hopDown(unsigned char from, unsigned char to) const;
        [[nodiscard]] SpinFermion transReduce() const;
        [[nodiscard]] inline int calcPeriod() const noexcept;
        inline void swap(SpinFermion& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] SpinlessFermion getSpinUp() const noexcept { return spinUp; }
        [[nodiscard]] SpinlessFermion getSpinDown() const noexcept { return spinDown; }
        [[nodiscard]] int getNumSite() const noexcept { return spinUp.getNumSite(); }
        [[nodiscard]] bool isVacuum() const noexcept { return spinUp.isVacuum() && spinDown.isVacuum(); }
        [[nodiscard]] bool isUpOccupy(unsigned char site) const noexcept { return spinUp.isOccupy(site); }
        [[nodiscard]] bool isDownOccupy(unsigned char site) const noexcept { return spinDown.isOccupy(site); }
        [[nodiscard]] unsigned int getNumSpinUpElectron() const noexcept { return spinUp.getNumElectron(); }
        [[nodiscard]] unsigned int getNumSpinDownElectron() const noexcept { return spinDown.getNumElectron(); }
        [[nodiscard]] inline unsigned int getNumPairedElectron() const noexcept;
    };

    inline SpinFermion::SpinFermion(SpinlessFermion spinUp_, SpinlessFermion spinDown_)
            : spinUp(spinUp_), spinDown(spinDown_) {}

    inline bool SpinFermion::operator>(const This& other) const noexcept {
        if (spinUp > other.spinUp)
            return true;
        if (spinUp == other.spinUp)
            return spinDown > other.spinDown;
        return false;
    }

    inline bool SpinFermion::operator<(const This& other) const noexcept {
        if (spinUp < other.spinUp)
            return true;
        if (spinUp == other.spinUp)
            return spinDown < other.spinDown;
        return false;
    }

    inline int SpinFermion::calcPeriod() const noexcept {
        const int result = lcm<int, false>(spinUp.calcPeriod(), spinDown.calcPeriod());
        assert(0 < result && result <= getNumSite() && "[Error]: Unexpected period, this is a bug");
        return result;
    }

    inline void SpinFermion::swap(SpinFermion& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        spinUp.swap(obj.spinUp);
        spinDown.swap(obj.spinDown);
    }

    inline unsigned int SpinFermion::getNumPairedElectron() const noexcept {
        return countOnes(spinUp.getOccupyBits() & spinDown.getOccupyBits());
    }

    inline std::ostream& operator<<(std::ostream& os, SpinFermion e) {
        return os << e.getSpinUp() << ' ' << e.getSpinDown();
    }
}
