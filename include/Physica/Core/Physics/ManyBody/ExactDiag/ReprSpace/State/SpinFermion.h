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
    template<unsigned int Dim, unsigned int NumSite>
    class SpinFermion : public State<SpinFermion<Dim, NumSite>> {
        using This = SpinFermion<Dim, NumSite>;
        using Base = State<This>;
        using SpinlessType = SpinlessFermion<Dim, NumSite>;
    private:
        SpinlessType spinUp;
        SpinlessType spinDown;
    public:
        SpinFermion() = default;
        inline SpinFermion(SpinlessType spinUp_, SpinlessType spinDown_);
        SpinFermion(const This&) = default;
        SpinFermion(This&&) noexcept = default;
        ~SpinFermion() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] bool operator==(const This& other) const noexcept { return (spinUp == other.spinUp) && (spinDown == other.spinDown); }
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
        [[nodiscard]] inline int hopUpSign(unsigned char from, unsigned char to) const;
        [[nodiscard]] inline int hopDownSign(unsigned char from, unsigned char to) const;
        [[nodiscard]] SpinFermion transReduce() const;
        [[nodiscard]] inline int lShiftSign() const;
        [[nodiscard]] inline int calcPeriod() const noexcept;
        inline void swap(SpinFermion& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] SpinlessType getSpinUp() const noexcept { return spinUp; }
        [[nodiscard]] SpinlessType getSpinDown() const noexcept { return spinDown; }
        [[nodiscard]] int getNumSite() const noexcept { return spinUp.getNumSite(); }
        [[nodiscard]] bool isVacuum() const noexcept { return spinUp.isVacuum() && spinDown.isVacuum(); }
        [[nodiscard]] bool isUpOccupy(unsigned char site) const noexcept { return spinUp.isOccupy(site); }
        [[nodiscard]] bool isDownOccupy(unsigned char site) const noexcept { return spinDown.isOccupy(site); }
        [[nodiscard]] unsigned int getNumParticle() const noexcept { return getNumSpinUpParticle() + getNumSpinDownParticle(); }
        [[nodiscard]] unsigned int getNumSpinUpParticle() const noexcept { return spinUp.getNumParticle(); }
        [[nodiscard]] unsigned int getNumSpinDownParticle() const noexcept { return spinDown.getNumParticle(); }
        [[nodiscard]] inline unsigned int getNumDoubleOccupy() const noexcept;
    };

    template<unsigned int Dim, unsigned int NumSite>
    inline std::ostream& operator<<(std::ostream& os, SpinFermion<Dim, NumSite> e) {
        return os << e.getSpinUp() << ' ' << e.getSpinDown();
    }
}

namespace Physica {
    template<unsigned int I1, unsigned int I2>
    class Traits<Core::SpinFermion<I1, I2>> {
    public:
        constexpr static unsigned int Dim = I1;
        constexpr static unsigned int NumSite = I2;
        constexpr static unsigned int SiteDOF = 4;
    };
}

#include "StateImpl/SpinFermionImpl.h"
