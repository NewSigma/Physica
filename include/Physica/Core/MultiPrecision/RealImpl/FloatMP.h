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

namespace Physica {
    template<>
    class Traits<Core::Real<Core::FloatMP>> {
    public:
        constexpr static Core::ScalarOption Option = Core::FloatMP;
        constexpr static int Order = 0;
        constexpr static bool isComplex = false;
        constexpr static bool isDifferentiable = false;
        constexpr static bool isForwardDiff = false;
        constexpr static bool isReverseDiff = false;

        using ValueType = Core::Real<Core::FloatMP>;
        using ScalarType = ValueType;
        using PtrTy = ScalarType*;
        using ConstPtrTy = const ScalarType*;
        using RefTy = ScalarType&;
        using ConstRefTy = const ScalarType&;
        using RealType = ScalarType;
        using ComplexType = Core::Complex<ScalarType>;
        using MachineType = ScalarType;
    };
}

namespace Physica::Core {
    constexpr int GlobalPrecision = 4;

    template<>
    class PHYSICA_API Real<FloatMP> : public ScalarBase<Real<FloatMP>> {
        static_assert(GlobalPrecision > 1, "GlobalPrecision must be larger than 1.");
    public:
        using ScalarType = Real<FloatMP>;
    private:
        //Store effective digits using little endian standard.
        MPUnit* __restrict byte;
        /**
         * Length of byte = abs(length).
         * sign of length and sign of Real are same. (when Real != 0)
         *
         * Warning: length can not equal to INT_MIN, or length will not return the correct answer.
         *
         * Optimize: use the end position of byte instead of length may improve performance.
         */
        int length;
        /**
         * Number = (x0 +- a * (2 ^ MPUnitWidth) ^ (1 - length)) * (2 ^ MPUnitWidth) ^ power
         *
         * FixIt: We have not considered overflow of power in our codes elsewhere.
         */
        int power;
    public:
        Real();
        Real(int length_, int power_);
        Real(int i);
        Real(SignedMPUnit unit);
        Real(double d);
        Real(const Integer& i);
        Real(const Rational& r);
        explicit Real(const char* s);
        explicit Real(const wchar_t* s);
        Real(const Real& s);
        Real(Real&& s) noexcept;
        ~Real();
        /* Operators */
        Real& operator=(Real s) noexcept { swap(s); return *this; }
        [[nodiscard]] MPUnit operator[](unsigned int index) const;
        [[nodiscard]] explicit operator double() const;
        [[nodiscard]] Real operator+(const Real& s) const;
        [[nodiscard]] Real operator-(const Real& s) const;
        [[nodiscard]] Real operator*(const Real& s) const;
        [[nodiscard]] Real operator/(const Real& s) const;
        [[nodiscard]] Real operator<<(int bits) const;
        [[nodiscard]] Real operator>>(int bits) const;
        [[nodiscard]] Real operator-() const;
        /* Operations */
        Real& toOpposite() noexcept { length = -length; return *this; }
        Real& toAbs() noexcept { length = getSize(); return *this; }
        void swap(Real& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] constexpr static ScalarOption getOption() { return FloatMP; }
        [[nodiscard]] int getLength() const noexcept { return length; }
        [[nodiscard]] int getPower() const noexcept { return power; }
        [[nodiscard]] int getSize() const noexcept { return std::abs(length); }
        [[nodiscard]] bool isZero() const { return byte[getSize() - 1] == 0; }
        [[nodiscard]] bool isPositive() const { return !isZero() && length > 0; }
        [[nodiscard]] bool isNegative() const { return !isZero() && length < 0; }
        [[nodiscard]] bool isInteger() const { return getSize() - 1 == power; }
        /* Setters */
        void setPower(int i) noexcept { power = i; }
        void setByte(unsigned int index, MPUnit value) { assert(index < static_cast<unsigned int>(getSize())); byte[index] = value; }
        /* Static members */
        static inline bool matchSign(const Real& s1, const Real& s2);
    private:
        /**
         * Degigned for performance,
         * this constructor should only be called by add(), sub() and etc.
         *
         * \param byte
         * byte must be allocated by malloc()
         */
        Real(MPUnit* byte_, int length_, int power_) : byte(byte_), length(length_), power(power_) {}
        /* Operations */
        inline void cutZero();
        /* Friends */
        friend class Integer;
        template<ScalarOption Option> __host__ __device__ friend Real<Option> square(const Real<Option>& s) noexcept;
        template<ScalarOption Option> __host__ __device__ friend Real<Option> sqrt(const Real<Option>& s) noexcept;
        template<ScalarOption Option> friend Real<Option> ln(const Real<Option>& s) noexcept;
        /* Static members */
        inline static Real<FloatMP> add(const Real& s1, const Real& s2);
        inline static Real<FloatMP> sub(const Real& s1, const Real& s2);
        inline static Real<FloatMP> mul(const Real& s1, const Real& s2);
        inline static Real<FloatMP> div(const Real& s1, const Real& s2);
        inline static bool cutLength(Real<FloatMP>& s);
    };

    inline Real<FloatMP>& operator++(Real<FloatMP>& s);
    inline Real<FloatMP>& operator--(Real<FloatMP>& s);
    inline Real<FloatMP> operator++(Real<FloatMP>& s, int);
    inline Real<FloatMP> operator--(Real<FloatMP>& s, int);
    /* Compare */
    bool absCompare(const Real<FloatMP>& s1, const Real<FloatMP>& s2);
    bool operator> (const Real<FloatMP>& s1, const Real<FloatMP>& s2);
    bool operator< (const Real<FloatMP>& s1, const Real<FloatMP>& s2);
    bool operator== (const Real<FloatMP>& s1, const Real<FloatMP>& s2);
}

namespace std {
    template<>
    struct numeric_limits<Physica::Core::Real<Physica::Core::FloatMP>> {
        using T = Physica::Core::Real<Physica::Core::FloatMP>;
    public:
        static T epsilon() noexcept {
            auto result = T(static_cast<Physica::Core::SignedMPUnit>(1));
            result.setPower(1 - Physica::Core::GlobalPrecision);
            return result;
        }
    };
}

#include "FloatMPImpl/Const.h"
#include "FloatMPImpl/FloatMPImpl.h"
#include "MathFPMP.h"
