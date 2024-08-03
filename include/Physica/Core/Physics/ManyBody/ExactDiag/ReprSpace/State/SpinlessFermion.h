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

#include <iostream>
#include "State.h"
#include "Physica/Core/Math/Discrete/Combination.h"

namespace Physica::Core {
    template<unsigned int Dim, unsigned int NumSite>
    class SpinlessFermion : public State<SpinlessFermion<Dim, NumSite>> {
        using This = SpinlessFermion<Dim, NumSite>;
        using Base = State<This>;
    public:
        using IntType = int64_t;
        static_assert(NumSite < sizeof(IntType) * CHAR_BIT, "[Error]: Unexpected large site number");
    private:
        IntType occupyBits;
    public:
        SpinlessFermion() = default;
        inline SpinlessFermion(IntType occupyBits_);
        SpinlessFermion(const This&) = default;
        SpinlessFermion(This&&) noexcept = default;
        ~SpinlessFermion() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] inline bool operator==(const This& other) const noexcept;
        [[nodiscard]] bool operator!=(const This& other) const noexcept { return !(*this == other); }
        [[nodiscard]] bool operator>(const This& other) const noexcept { return occupyBits > other.occupyBits; }
        [[nodiscard]] bool operator<(const This& other) const noexcept { return occupyBits < other.occupyBits; }
        [[nodiscard]] bool operator>=(const This& other) const noexcept { return occupyBits >= other.occupyBits; }
        [[nodiscard]] bool operator<=(const This& other) const noexcept { return occupyBits <= other.occupyBits; }
        [[nodiscard]] inline This operator<<(int shift) const noexcept;
        [[nodiscard]] inline This operator>>(int shift) const noexcept;
        void operator<<=(int shift) noexcept { (*this) = (*this) << shift; }
        void operator>>=(int shift) noexcept { (*this) = (*this) >> shift; }
        /* Operations */
        [[nodiscard]] inline SpinlessFermion hop(uint8_t from, uint8_t to) const;
        [[nodiscard]] inline int hopSign(uint8_t from, uint8_t to) const;
        [[nodiscard]] SpinlessFermion transReduce(int period = 1) const;
        [[nodiscard]] int lShiftSign() const;
        [[nodiscard]] int calcPeriod() const;
        inline void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] IntType getOccupyBits() const noexcept { return occupyBits; }
        [[nodiscard]] bool isVacuum() const noexcept { return occupyBits == 0; }
        [[nodiscard]] inline bool isOccupy(uint8_t site) const noexcept;
        [[nodiscard]] unsigned int getNumParticle() const noexcept { return countOnes(occupyBits); }
        [[nodiscard]] inline bool isTransReducible(int period = 1) const noexcept;
        /* Static members */
        template<class RandomGenerator>
        [[nodiscard]] static This random_state(RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] static This random_state(size_t numParticle, RandomGenerator& gen);
        [[nodiscard]] constexpr static IntType makeFullMask() noexcept { return (static_cast<IntType>(1) << NumSite) - 1; }
        [[nodiscard]] constexpr static IntType makeHighMask() noexcept { return (static_cast<IntType>(1) << (NumSite - 1)); }
    };

    template<unsigned int Dim, unsigned int NumSite>
    std::ostream& operator<<(std::ostream& os, SpinlessFermion<Dim, NumSite> e) {
        auto mask = e.makeHighMask();
        for (int i = 0; i < int(NumSite); ++i) {
            const bool flag = (e.getOccupyBits() & mask) == 0;
            os << (flag ? '0' : '1');
            mask >>= 1;
        }
        return os;
    }
}

namespace Physica {
    template<unsigned int I1, unsigned int I2>
    class Traits<Core::SpinlessFermion<I1, I2>> {
    public:
        constexpr static unsigned int Dim = I1;
        constexpr static unsigned int NumSite = I2;
        constexpr static unsigned int SiteDOF = 2;

        static_assert(1 <= Dim && Dim <= 3, "[Error]: Invalid Dim");
        static_assert(0 < NumSite && NumSite < 64, "[Error]: Invalid site number");
    };
}

namespace std {
    template<unsigned int Dim, unsigned int NumSite>
    struct hash<Physica::Core::SpinlessFermion<Dim, NumSite>> {
        using T = Physica::Core::SpinlessFermion<Dim, NumSite>;

        std::size_t operator()(const T& psi) const noexcept {
            return std::hash<typename T::IntType>{}(psi.getOccupyBits());
        }
        /**
         * Reference:
         * [1] Comput. Phys. Commun. 92(1), 11-15 (1995); https://doi.org/10.1016/0010-4655(95)00108-R
         * [2] Comput. Phys. Commun. 224, 81-89 (2018); https://doi.org/10.1016/j.cpc.2017.11.011
         */
        std::size_t perfect_hash(const T& psi) const noexcept {
            auto bits = psi.getOccupyBits();
            int i = 0;
            size_t result = 0;
            for (int site = 0; site < int(NumSite); ++site) {
                const bool isOccupy = bits % 2 == 1;
                if (isOccupy) {
                    i += 1;
                    if (site >= i)
                        result += Physica::Core::combination(site, i);
                }
                bits >>= 1;
            }
            return result;
        }
    };
}

#include "StateImpl/SpinlessFermionImpl.h"
