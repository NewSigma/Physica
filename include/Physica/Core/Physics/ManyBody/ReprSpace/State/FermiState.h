/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "SpinState.h"

namespace Physica {
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
        FermiState(Spin spinUp_, Spin spinDown_);
        FermiState(const This&) = default;
        FermiState(This&&) noexcept = default;
        ~FermiState() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] bool operator==(const This& other) const noexcept;
        [[nodiscard]] auto operator<=>(const This& other) const noexcept;
        [[nodiscard]] This operator<<(int shift) const noexcept { return This(spinUp << shift, spinDown << shift); }
        [[nodiscard]] This operator>>(int shift) const noexcept { return This(spinUp >> shift, spinDown >> shift); }
        void operator<<=(int shift) noexcept { (*this) = (*this) << shift; }
        void operator>>=(int shift) noexcept { (*this) = (*this) >> shift; }
        /* Operations */
        [[nodiscard]] This hopUp(uint8_t from, uint8_t to) const noexcept;
        [[nodiscard]] This hopDown(uint8_t from, uint8_t to) const noexcept;
        [[nodiscard]] int hopUpSign(uint8_t from, uint8_t to) const noexcept;
        [[nodiscard]] int hopDownSign(uint8_t from, uint8_t to) const noexcept;
        [[nodiscard]] This transReduce() const;
        [[nodiscard]] int lShiftSign() const;
        [[nodiscard]] int calcPeriod() const noexcept;
        void swap(This& __restrict obj) noexcept;
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
        [[nodiscard]] int getNumDoubleOccupy() const noexcept;
        /* Static member */
        template<RNG R>
        [[nodiscard]] static This random_state();
        template<RNG R>
        [[nodiscard]] static This random_state(size_t numSpinUp, size_t numSpinDown);
    };

    template<int Dim, int NumSite>
    std::ostream& operator<<(std::ostream& os, FermiState<Dim, NumSite> e) {
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
