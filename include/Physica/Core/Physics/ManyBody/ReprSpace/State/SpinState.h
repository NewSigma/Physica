/*
 * Copyright 2024-2026 Weibo He.
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

#include "Physica/Core/Math/Discrete/Combination.h"
#include "Physica/Core/Math/Random/Random.h"
#include "State.h"

namespace Physica {
    template<int Dim, int NumSite>
    class SpinState : public State<SpinState<Dim, NumSite>> {
        using This = SpinState<Dim, NumSite>;
        using Base = State<This>;
    public:
        using IntType = uint64_t;
        static_assert(NumSite <= sizeof(IntType) * CHAR_BIT, "[Error]: Unexpected large site number");
    private:
        IntType occupyBits = 0;
    public:
        SpinState() = default;
        SpinState(IntType occupyBits_) noexcept;
        SpinState(const This&) = default;
        SpinState(This&&) noexcept = default;
        ~SpinState() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] bool operator==(const This& other) const noexcept;
        [[nodiscard]] auto operator<=>(const This& other) const noexcept;
        [[nodiscard]] bool operator[](int8_t site) const noexcept { return isOccupy(site); }
        [[nodiscard]] This operator<<(int shift) const noexcept;
        [[nodiscard]] This operator>>(int shift) const noexcept;
        [[nodiscard]] This operator&(const This& psi) const noexcept;
        [[nodiscard]] This operator|(const This& psi) const noexcept;
        [[nodiscard]] This operator^(const This& psi) const noexcept;
        [[nodiscard]] This operator~() const noexcept;
        void operator<<=(int shift) noexcept { (*this) = (*this) << shift; }
        void operator>>=(int shift) noexcept { (*this) = (*this) >> shift; }
        void operator&=(const This& psi) noexcept { (*this) = (*this) & psi; }
        void operator|=(const This& psi) noexcept { (*this) = (*this) | psi; }
        void operator^=(const This& psi) noexcept { (*this) = (*this) ^ psi; }
        /* Operations */
        [[nodiscard]] This flip(int8_t site) const;
        [[nodiscard]] This hop(int8_t from, int8_t to) const noexcept;
        [[nodiscard]] int hopSign(int8_t from, int8_t to) const noexcept;
        [[nodiscard]] This transReduce(int period = 1) const;
        [[nodiscard]] int lShiftSign() const;
        [[nodiscard]] int calcPeriod() const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] IntType getOccupyBits() const noexcept { return occupyBits; }
        [[nodiscard]] bool isVacuum() const noexcept { return occupyBits == 0; }
        [[nodiscard]] bool isOccupy(int8_t site) const noexcept;
        [[nodiscard]] int getNumParticle() const noexcept { return std::popcount(occupyBits); }
        [[nodiscard]] bool isTransReducible(int period = 1) const noexcept;

        [[nodiscard]] bool isSpinUp(int8_t site) const noexcept;
        [[nodiscard]] bool isSpinDown(int8_t site) const noexcept;
        [[nodiscard]] int getNumSpinUp() const noexcept;
        [[nodiscard]] int getNumSpinDown() const noexcept;
        /* Static members */
        template<RNG R>
        [[nodiscard]] static This random_state();
        template<RNG R>
        [[nodiscard]] static This random_state(size_t numParticle);
        [[nodiscard]] constexpr static IntType makeFullMask() noexcept { return (static_cast<IntType>(1) << NumSite) - 1; }
        [[nodiscard]] constexpr static IntType makeHighMask() noexcept { return (static_cast<IntType>(1) << (NumSite - 1)); }
    };

    template<int Dim, int NumSite>
    std::ostream& operator<<(std::ostream& os, SpinState<Dim, NumSite> e) {
        auto mask = e.makeHighMask();
        for (int i = 0; i < NumSite; ++i) {
            const bool flag = (e.getOccupyBits() & mask) == 0;
            os << (flag ? '0' : '1');
            mask >>= 1;
        }
        return os;
    }
}

namespace Physica {
    template<int I1, int I2>
    class Traits<SpinState<I1, I2>> {
    public:
        constexpr static int Dim = I1;
        constexpr static int NumSite = I2;
        constexpr static int SiteDOF = 2;

        static_assert(NumSite < 64, "[Error]: Use Dynamic if NumSite is large");
    };
}

namespace std {
    template<int Dim, int NumSite>
    struct hash<Physica::SpinState<Dim, NumSite>> {
        using T = Physica::SpinState<Dim, NumSite>;

        std::size_t operator()(const T& psi) const noexcept {
            return std::hash<typename T::IntType>{}(psi.getOccupyBits());
        }
        /**
         * Reference:
         * [1] Comput. Phys. Commun. 92(1), 11-15 (1995); https://doi.org/10.1016/0010-4655(95)00108-R
         * [2] Comput. Phys. Commun. 224, 81-89 (2018); https://doi.org/10.1016/j.cpc.2017.11.011
         */
        static std::size_t perfect_hash(const T& psi) noexcept {
            auto bits = psi.getOccupyBits();
            int i = 0;
            size_t result = 0;
            for (int site = 0; site < NumSite; ++site) {
                const bool isOccupy = bits % 2 == 1;
                if (isOccupy) {
                    i += 1;
                    if (site >= i)
                        result += Physica::binomial(site, i);
                }
                bits >>= 1;
            }
            return result;
        }
    };
}

#include "StateImpl/SpinStateImpl.h"
