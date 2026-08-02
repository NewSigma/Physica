/*
 * Copyright 2020-2026 Weibo He.
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

#include <cmath>
#include <cassert>
#include "Physica/Core/Utils/Allocator/HostAllocator.h"
#include "Scalar.h"

namespace Physica {
    template<FloatPrec Prec> class Real;

    class PHYSICA_API Integer {
        using This = Integer;
        //Store effective digits using little endian standard.
        MPUnit* __restrict byte = nullptr;
        /*
         * Length of byte = abs(length).
         * sign of length and sign of Integer are same. (when Integer != 0)
         *
         * Warning: length can not equal to INT_MIN, or length will not return the correct answer.
         *
         * Optimize: use the end position of byte instead of length may improve performance.
         * Optimize: length of byte is known and is stored in stack so the extra space that is taken up by malloc() is unnecessary.
         */
        int length;
    public:
        Integer(int i); //NOLINT Conversion is always available.
        template<FloatPrec Prec>
        explicit Integer(const Real<Prec>& s);
        Integer(const This& obj);
        Integer(This&& obj) noexcept;
        ~Integer();
        /* Operators */
        Integer& operator=(Integer obj) noexcept;
        Integer operator+(const Integer& i) const;
        Integer operator-(const Integer& i) const;
        Integer operator*(const Integer& i) const;
        Integer operator/(const Integer& i) const;
        Integer operator%(const Integer& i) const;
        Integer operator<<(int bits) const;
        Integer operator>>(int bits) const;
        Integer operator-() const;
        bool operator>(const Integer& i) const;
        bool operator<(const Integer& i) const;
        bool operator==(const Integer& i) const;
        bool operator!=(const Integer& i) const { return !(*this == i); }
        bool operator>=(const Integer& i) const { return !(*this < i); }
        bool operator<=(const Integer& i) const { return !(*this > i); }
        void operator+=(const Integer& i) { *this = *this + i; }
        void operator-=(const Integer& i) { *this = *this - i; }
        void operator*=(const Integer& i) { *this = *this * i; }
        void operator/=(const Integer& i) { *this = *this / i; }
        void operator%=(const Integer& i) { *this = *this % i; }
        void operator<<=(int bits) { *this = *this << bits; }
        void operator>>=(int bits) { *this = *this >> bits; }
        Integer& operator++();
        Integer& operator--();
        Integer operator++(int);
        Integer operator--(int);
        explicit operator double() const;
        /* Operations */
        Integer& toOpposite() noexcept { length = -length; return *this; }
        Integer& toAbs() noexcept { length = getSize(); return *this; }
        void swap(Integer& __restrict i) noexcept;
        /* Getters */
        [[nodiscard]] const MPUnit* getByte() const noexcept { return byte; }
        [[nodiscard]] int getLength() const noexcept { return length; }
        [[nodiscard]] int getSize() const noexcept { return std::abs(length); }
        [[nodiscard]] bool isZero() const { return byte[getSize() - 1] == 0; }
        [[nodiscard]] bool isPositive() const { return !isZero() && length > 0; }
        [[nodiscard]] bool isNegative() const { return !isZero() && length < 0; }
        [[nodiscard]] bool isOdd() const { return byte[0] & 1U; }
        [[nodiscard]] bool isEven() const { return !(byte[0] & 1U); }
        /* Setters */
        void setSign(bool sign) noexcept { length = sign ? length : -length; }
        void setByte(unsigned int index, MPUnit value) { assert(index < static_cast<unsigned int>(getSize())); byte[index] = value; }
        /* Static members */
        static inline bool matchSign(const Integer& i1, const Integer& i2);
        static bool absCompare(const Integer& i1, const Integer& i2);
    protected:
        /**
         * Designed for performance,
         * this constructor should only be called operator+, -, *, / and etc.
         */
        Integer(MPUnit* byte_, int length_) : byte(byte_), length(length_) {}
        /* Operations */
        void cutZero();
        /* Static members */
        static Integer integerAddImpl(const Integer& i1, const Integer& i2);
        static Integer integerSubImpl(const Integer& i1, const Integer& i2);
    };

    template<FloatPrec Prec>
    Integer::Integer(const Real<Prec>& s) {
        if constexpr (Prec == FloatMP) {
            const auto power = s.getPower();
            if (power < 0) {
                byte = HostAllocator<MPUnit>{}.allocate(1);
                byte[0] = 0;
                length = 1;
            }
            length = power + 1;
            byte = HostAllocator<MPUnit>{}.allocate(length);
            memcpy(byte, s.byte, length * sizeof(MPUnit));
        }
        else
            *this = Integer(static_cast<int>(s.toMachine()));
    }
    /**
     * Returns true if i1 and i2 has the same sign. Both i1 and i2 do not equal to zero.
     * This function provide a quick sign check compare to using isPositive() and isNegative().
     */
    inline bool Integer::matchSign(const Integer& i1, const Integer& i2) {
        assert(!i1.isZero() && !i2.isZero());
        return (i1.length ^ i2.length) >= 0;
    }

    Integer factorial(const Integer& i);

    inline void swap(Integer& i1, Integer& i2) noexcept {
        i1.swap(i2);
    }
}

namespace std {
    template<>
    struct numeric_limits<Physica::Integer> {
        constexpr static bool is_integer = true;
    };
}
