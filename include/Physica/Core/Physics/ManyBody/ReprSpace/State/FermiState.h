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

#include "SpinState.h"

namespace Physica::Core {
    template<int Dim, int NumSite>
    class FermiState : public State<FermiState<Dim, NumSite>> {
        using This = FermiState<Dim, NumSite>;
        using Base = State<This>;
        using Spin = SpinState<Dim, NumSite>;
    public:
        using IntType = Spin::IntType;
    private:
        Spin spinUp;
        Spin spinDown;
    public:
        FermiState() = default;
        inline FermiState(Spin spinUp_, Spin spinDown_);
        FermiState(const This&) = default;
        FermiState(This&&) noexcept = default;
        ~FermiState() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] bool operator==(const This& other) const noexcept { return (spinUp == other.spinUp) && (spinDown == other.spinDown); }
        [[nodiscard]] bool operator!=(const This& other) const noexcept { return !(*this == other); }
        [[nodiscard]] inline bool operator>(const This& other) const noexcept;
        [[nodiscard]] inline bool operator<(const This& other) const noexcept;
        [[nodiscard]] bool operator>=(const This& other) const noexcept { return !((*this) < other); }
        [[nodiscard]] bool operator<=(const This& other) const noexcept { return !((*this) > other); }
        [[nodiscard]] This operator<<(int shift) const noexcept { return FermiState(spinUp << shift, spinDown << shift); }
        [[nodiscard]] This operator>>(int shift) const noexcept { return FermiState(spinUp >> shift, spinDown >> shift); }
        void operator<<=(int shift) noexcept { (*this) = (*this) << shift; }
        void operator>>=(int shift) noexcept { (*this) = (*this) >> shift; }
        /* Operations */
        [[nodiscard]] FermiState hopUp(uint8_t from, uint8_t to) const;
        [[nodiscard]] FermiState hopDown(uint8_t from, uint8_t to) const;
        [[nodiscard]] inline int hopUpSign(uint8_t from, uint8_t to) const;
        [[nodiscard]] inline int hopDownSign(uint8_t from, uint8_t to) const;
        [[nodiscard]] FermiState transReduce() const;
        [[nodiscard]] inline int lShiftSign() const;
        [[nodiscard]] inline int calcPeriod() const noexcept;
        inline void swap(FermiState& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] Spin getSpinUp() const noexcept { return spinUp; }
        [[nodiscard]] Spin getSpinDown() const noexcept { return spinDown; }
        [[nodiscard]] int getNumSite() const noexcept { return spinUp.getNumSite(); }
        [[nodiscard]] bool isVacuum() const noexcept { return spinUp.isVacuum() && spinDown.isVacuum(); }
        [[nodiscard]] bool isUpOccupy(uint8_t site) const noexcept { return spinUp.isOccupy(site); }
        [[nodiscard]] bool isDownOccupy(uint8_t site) const noexcept { return spinDown.isOccupy(site); }
        [[nodiscard]] int getNumParticle() const noexcept { return getNumSpinUpParticle() + getNumSpinDownParticle(); }
        [[nodiscard]] int getNumSpinUpParticle() const noexcept { return spinUp.getNumParticle(); }
        [[nodiscard]] int getNumSpinDownParticle() const noexcept { return spinDown.getNumParticle(); }
        [[nodiscard]] inline int getNumDoubleOccupy() const noexcept;
        /* Static member */
        template<RandomGenerator R>
        [[nodiscard]] static This random_state();
        template<RandomGenerator R>
        [[nodiscard]] static This random_state(size_t numSpinUp, size_t numSpinDown);
    };

    template<int Dim, int NumSite>
    inline std::ostream& operator<<(std::ostream& os, FermiState<Dim, NumSite> e) {
        return os << e.getSpinUp() << ' ' << e.getSpinDown();
    }
}

namespace Physica {
    template<int I1, int I2>
    class Traits<FermiState<I1, I2>> {
    public:
        constexpr static int Dim = I1;
        constexpr static int NumSite = I2;
        constexpr static int SiteDOF = 4;
    };
}

#include "StateImpl/FermiStateImpl.h"
