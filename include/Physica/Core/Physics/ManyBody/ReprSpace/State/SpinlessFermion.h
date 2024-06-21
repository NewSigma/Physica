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

#include <iostream>
#include "State.h"

namespace Physica::Core {
    template<unsigned int Dim, unsigned int NumSite> class SpinlessFermion;

    namespace Internal {
        template<unsigned int I1, unsigned int I2>
        class Traits<SpinlessFermion<I1, I2>> {
        public:
            constexpr static unsigned int Dim = I1;
            constexpr static unsigned int NumSite = I2;
            constexpr static unsigned int SiteDOF = 2;
        };
    }

    template<unsigned int Dim, unsigned int NumSite>
    class SpinlessFermion : public State<SpinlessFermion<Dim, NumSite>> {
        static_assert(1 <= Dim && Dim <= 3, "[Error]: Invalid Dim");
        static_assert(NumSite > 0, "[Error]: Invalid site number");
        using This = SpinlessFermion<Dim, NumSite>;
        using Base = State<This>;
    public:
        using typename Base::IndexType;
    private:
        uint64_t occupyBits;

        static_assert(NumSite < sizeof(occupyBits) * CHAR_BIT, "[Error]: Unexpected large site number");
    public:
        SpinlessFermion() = default;
        inline SpinlessFermion(uint64_t occupyBits_);
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
        [[nodiscard]] inline SpinlessFermion hop(unsigned char from, unsigned char to) const;
        [[nodiscard]] inline SpinlessFermion hop(IndexType dims, IndexType from, IndexType to) const;
        [[nodiscard]] int hopSign(unsigned char from, unsigned char to) const;
        [[nodiscard]] int hopSign(IndexType dims, IndexType from, IndexType to) const;
        [[nodiscard]] SpinlessFermion transReduce(int period = 1) const;
        [[nodiscard]] int calcPeriod() const;
        inline void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] uint64_t getOccupyBits() const noexcept { return occupyBits; }
        [[nodiscard]] bool isVacuum() const noexcept { return occupyBits == 0; }
        [[nodiscard]] inline bool isOccupy(unsigned char site) const noexcept;
        [[nodiscard]] unsigned int getNumElectron() const noexcept { return countOnes(occupyBits); }
        [[nodiscard]] inline bool isTransReducible(int period = 1) const noexcept;
        [[nodiscard]] inline uint64_t makeFullMask() const noexcept;
        [[nodiscard]] inline uint64_t makeHighMask() const noexcept;
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

#include "StateImpl/SpinlessFermionImpl.h"
