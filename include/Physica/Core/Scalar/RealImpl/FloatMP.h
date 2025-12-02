/*
 * Copyright 2024-2025 Weibo He.
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

#include "../Real.h"

namespace Physica {
    template<>
    class Traits<Real<FloatMP>> {
    public:
        constexpr static FloatPrec Prec = FloatMP;
        constexpr static int Order = 0;
        constexpr static bool isComplex = false;
        constexpr static bool isForwardDiff = false;
        constexpr static bool isReverseDiff = false;

        using ValueType = Real<FloatMP>;
        using ScalarType = ValueType;
        using PtrTy = ScalarType*;
        using ConstPtrTy = const ScalarType*;
        using RefTy = ScalarType&;
        using ConstRefTy = const ScalarType&;
        using RealType = ScalarType;
        using ComplexType = Complex<ScalarType>;
        using MachineType = double;
    };
}

namespace Physica {
    constexpr int GlobalPrecision = 4;
    static_assert(GlobalPrecision > 1, "GlobalPrecision must be larger than 1.");

    template<>
    class PHYSICA_API Real<FloatMP> : public ScalarBase<Real<FloatMP>>, public CRCoro<Real<FloatMP>> {
        using This = Real<FloatMP>;
        using Base = ScalarBase<This>;
    private:
        //Store effective digits using little endian standard.
        MPUnit* __restrict byte;
        /**
         * Length of byte = abs(length).
         * sign of length and sign of Real are same. (when Real != 0)
         *
         * Warning: length can not equal to INT_MIN, or it will lead the incorrect result.
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
        Real(std::initializer_list<MPUnit> bytes_, int length_, int power_);
        Real(SignedMPUnit x);
        Real(std::integral auto x);
        Real(double d) noexcept;
        Real(const Integer& i);
        Real(const Rational& r);
        template<Scalar T>
        Real(const T& x) requires(!Diffable<T>);
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
        [[nodiscard]] bool operator> (const Real& other) const;
        [[nodiscard]] bool operator< (const Real& other) const;
        [[nodiscard]] bool operator== (const Real& other) const;
        /* Operations */
        using Base::random_uniform;
        using Base::random_normal;
        Real& toOpposite() noexcept { length = -length; return *this; }
        Real& toAbs() noexcept { length = getSize(); return *this; }

        void swap(Real& __restrict obj) noexcept;
        void dump() const noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ int getLength() const noexcept { return length; }
        [[nodiscard]] __host__ __device__ int getPower() const noexcept { return power; }
        [[nodiscard]] __host__ __device__ int getSize() const noexcept { return std::abs(length); }
        [[nodiscard]] MachineType toMachine() const noexcept { return double(*this); }
        [[nodiscard]] bool isZero() const { return byte[getSize() - 1] == 0; }
        [[nodiscard]] constexpr static bool isSubNormal() noexcept { return false; }
        [[nodiscard]] bool isPositive() const { return !isZero() && length > 0; }
        [[nodiscard]] bool isNegative() const { return !isZero() && length < 0; }
        [[nodiscard]] bool isInteger() const { return getSize() - 1 == power; }
        /* Setters */
        void setPower(int i) noexcept { power = i; }
        void setByte(unsigned int index, MPUnit value) noexcept;
        /* Static members */
        template<RNG R>
        [[nodiscard]] static Real random_uniform();
        template<RNG R>
        [[nodiscard]] static Real random_normal();
        [[nodiscard]] static bool matchSign(const Real& s1, const Real& s2) noexcept;
        [[nodiscard]] static double toDouble(MPUnit* __restrict byte, int length, int power) noexcept;
    private:
        /**
         * Degigned for performance,
         * this constructor should only be called by add(), sub() and etc.
         *
         * \param byte
         * byte must be allocated by new
         */
        Real(MPUnit* byte_, int length_, int power_) : byte(byte_), length(length_), power(power_) {}
        /* Operations */
        void cutZero();
        /* Friends */
        friend class Integer;
        template<FloatPrec Prec> __host__ __device__ friend Real<Prec> square(const Real<Prec>& x) noexcept;
        template<FloatPrec Prec> __host__ __device__ friend Real<Prec> sqrt(const Real<Prec>& x) noexcept;
        template<FloatPrec Prec> __host__ __device__ friend Real<Prec> ln(const Real<Prec>& x) noexcept;
        /* Static members */
        static Real<FloatMP> add(const Real& s1, const Real& s2) noexcept;
        static Real<FloatMP> sub(const Real& s1, const Real& s2) noexcept;
        static Real<FloatMP> mul(const Real& s1, const Real& s2) noexcept;
        static Real<FloatMP> div(const Real& s1, const Real& s2) noexcept;
        static bool cutLength(Real<FloatMP>& s);
    };

    Real<FloatMP>::Real(std::integral auto x) : Real(SignedMPUnit(x)) {
        assert(x <= std::numeric_limits<SignedMPUnit>::max());
    }

    template<Scalar T>
    Real<FloatMP>::Real(const T& x) requires(!Diffable<T>) : Real(double(x)) {}

    template<RNG R>
    Real<FloatMP> Real<FloatMP>::random_uniform() {
        return Real<Float64>::random_uniform<R>();
    }

    template<RNG R>
    Real<FloatMP> Real<FloatMP>::random_normal() {
        return Real<Float64>::random_normal<R>();
    }

    inline void Real<FloatMP>::setByte(unsigned int index, MPUnit value) noexcept {
        assert(index < static_cast<unsigned int>(getSize()));
        byte[index] = value;
    }

    inline std::ostream& operator<<(std::ostream& os, const Real<FloatMP>& x) {
        return os << std::format("{}", x.toMachine());
    }

    PHYSICA_API Real<FloatMP>& operator++(Real<FloatMP>& s);
    PHYSICA_API Real<FloatMP>& operator--(Real<FloatMP>& s);
    PHYSICA_API Real<FloatMP> operator++(Real<FloatMP>& s, int);
    PHYSICA_API Real<FloatMP> operator--(Real<FloatMP>& s, int);
    PHYSICA_API bool absCompare(const Real<FloatMP>& s1, const Real<FloatMP>& s2);
}

#include "FloatMPImpl/Const.h"

namespace std {
    template<>
    struct numeric_limits<Physica::Real<Physica::FloatMP>> : public numeric_limits<double> {};

    template<>
    struct PHYSICA_API formatter<Physica::Real<Physica::FloatMP>, char> {
        constexpr static auto parse(std::format_parse_context& ctx) { return ctx.begin(); }
        static auto format(const Physica::Real<Physica::FloatMP>& obj, std::format_context& ctx) -> std::format_context::iterator;
    };
}
