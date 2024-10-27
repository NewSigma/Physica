/*
 * Copyright 2020-2023 Weibo He.
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
#include <utility>
#include "MultiPrecisionType.h"

namespace Physica::Core {
    class PHYSICA_API Integer {
        //Store effective digits using little endian standard.
        MPUnit* __restrict byte;
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
        template<ScalarOption Option>
        explicit Integer(const Scalar<Option>& s);
        explicit Integer(const Scalar<FloatMP>& s);
        Integer(const Integer& toCopy);
        Integer(Integer&& toMove) noexcept;
        ~Integer();
        /* Operators */
        Integer& operator=(const Integer& toCopy);
        Integer& operator=(Integer&& toMove) noexcept;
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
         * Degigned for performance,
         * this constructor should only be called operator+, -, *, / and etc.
         *
         * \param byte
         * byte must be allocated by malloc()
         */
        Integer(MPUnit* byte_, int length_) : byte(byte_), length(length_) {}
        /* Operations */
        void cutZero();
        /* Static members */
        static inline Integer integerAddImpl(const Integer& i1, const Integer& i2);
        static inline Integer integerSubImpl(const Integer& i1, const Integer& i2);
    };
}

namespace std {
    template<>
    struct numeric_limits<Physica::Core::Integer> {
        constexpr static bool is_integer = true;
    };

    template<>
    inline void swap(Physica::Core::Integer& __restrict i1, Physica::Core::Integer& __restrict i2) noexcept {
        i1.swap(i2);
    }
}

#include "IntegerImpl/IntegerImpl.h"
