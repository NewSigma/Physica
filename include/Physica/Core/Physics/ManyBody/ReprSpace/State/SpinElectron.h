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

#include "SpinlessElectron.h"
#include "Physica/Core/Math/NumberTheory/NumberTheory.h"

namespace Physica::Core {
    class SpinElectron;

    namespace Internal {
        template<>
        class Traits<SpinElectron> {
        public:
            constexpr static unsigned int SiteDOF = 4;
        };
    }

    class PHYSICA_API SpinElectron : public State<SpinElectron> {
        using This = SpinElectron;

        SpinlessElectron spinUp;
        SpinlessElectron spinDown;
    public:
        SpinElectron() = default;
        inline SpinElectron(SpinlessElectron spinUp_, SpinlessElectron spinDown_);
        SpinElectron(const SpinElectron&) = default;
        SpinElectron(SpinElectron&&) noexcept = default;
        ~SpinElectron() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] bool operator==(const This& other) const noexcept { return spinUp == other.spinUp && spinDown == other.spinDown; }
        [[nodiscard]] bool operator!=(const This& other) const noexcept { return !(*this == other); }
        [[nodiscard]] inline bool operator>(const This& other) const noexcept;
        [[nodiscard]] inline bool operator<(const This& other) const noexcept;
        [[nodiscard]] bool operator>=(const This& other) const noexcept { return !((*this) < other); }
        [[nodiscard]] bool operator<=(const This& other) const noexcept { return !((*this) > other); }
        [[nodiscard]] This operator<<(int shift) const noexcept { return SpinElectron(spinUp << shift, spinDown << shift); }
        [[nodiscard]] This operator>>(int shift) const noexcept { return SpinElectron(spinUp >> shift, spinDown >> shift); }
        void operator<<=(int shift) noexcept { (*this) = (*this) << shift; }
        void operator>>=(int shift) noexcept { (*this) = (*this) >> shift; }
        /* Operations */
        [[nodiscard]] SpinElectron hopUp(unsigned char from, unsigned char to) const;
        [[nodiscard]] SpinElectron hopDown(unsigned char from, unsigned char to) const;
        [[nodiscard]] SpinElectron transReduce() const;
        [[nodiscard]] inline int calcPeriod() const noexcept;
        inline void swap(SpinElectron& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] SpinlessElectron getSpinUp() const noexcept { return spinUp; }
        [[nodiscard]] SpinlessElectron getSpinDown() const noexcept { return spinDown; }
        [[nodiscard]] int getNumSite() const noexcept { return spinUp.getNumSite(); }
        [[nodiscard]] bool isVacuum() const noexcept { return spinUp.isVacuum() && spinDown.isVacuum(); }
        [[nodiscard]] bool isUpOccupy(unsigned char site) const noexcept { return spinUp.isOccupy(site); }
        [[nodiscard]] bool isDownOccupy(unsigned char site) const noexcept { return spinDown.isOccupy(site); }
        [[nodiscard]] unsigned int getNumSpinUpElectron() const noexcept { return spinUp.getNumElectron(); }
        [[nodiscard]] unsigned int getNumSpinDownElectron() const noexcept { return spinDown.getNumElectron(); }
        [[nodiscard]] inline unsigned int getNumPairedElectron() const noexcept;
    };

    inline SpinElectron::SpinElectron(SpinlessElectron spinUp_, SpinlessElectron spinDown_)
            : spinUp(spinUp_), spinDown(spinDown_) {}

    inline bool SpinElectron::operator>(const This& other) const noexcept {
        if (spinUp > other.spinUp)
            return true;
        if (spinUp == other.spinUp)
            return spinDown > other.spinDown;
        return false;
    }

    inline bool SpinElectron::operator<(const This& other) const noexcept {
        if (spinUp < other.spinUp)
            return true;
        if (spinUp == other.spinUp)
            return spinDown < other.spinDown;
        return false;
    }

    inline int SpinElectron::calcPeriod() const noexcept {
        const int result = lcm<int, false>(spinUp.calcPeriod(), spinDown.calcPeriod());
        assert(0 < result && result <= getNumSite() && "[Error]: Unexpected period, this is a bug");
        return result;
    }

    inline void SpinElectron::swap(SpinElectron& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        spinUp.swap(obj.spinUp);
        spinDown.swap(obj.spinDown);
    }

    inline unsigned int SpinElectron::getNumPairedElectron() const noexcept {
        return countOnes(spinUp.getOccupyBits() & spinDown.getOccupyBits());
    }

    inline std::ostream& operator<<(std::ostream& os, SpinElectron e) {
        return os << e.getSpinUp() << ' ' << e.getSpinDown();
    }
}
