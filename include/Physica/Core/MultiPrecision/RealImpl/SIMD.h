/*
 * Copyright 2023-2024 Weibo He.
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

#include <vectorclass/vectorclass.h>
#include <vectorclass/vectormath_exp.h>
#include "Physica/PlainStruct.h"
#include "Physica/Core/Utils/Container/Array.h"
#include "SIMDImpl/Instruset.h"

namespace Physica::Core {
    template<Scalar T, size_t Size> class BoolSIMD;

    template<Scalar T, size_t Size>
    class SIMD : private Traits<SIMD<T, Size>>::BaseType {
        using This = SIMD<T, Size>;
        using Base = Traits<This>::BaseType;
        using ValueType = T::ValueType;
        using MachineType = T::MachineType;
        using HalfType = std::conditional<sizeof(Base) * CHAR_BIT != 128, SIMD<T, Size / 2>, PlainStruct<void>>::type;
        constexpr static bool isForward = T::isForwardDiff;
    public:
        using BoolSIMDType = Traits<This>::BoolSIMDType;
        using PlainPacket = This;

        constexpr static bool isSeparatable = !std::is_same<HalfType, PlainStruct<void>>::value;
    public:
        SIMD() = default;
        explicit SIMD(T s) : Base(s.toMachine()) {}
        SIMD(T s, int count);
        SIMD(Base value) : Base(value) {}
        SIMD(HalfType a, HalfType b);
        using Base::Base;
        SIMD(const SIMD&) = default;
        SIMD(SIMD&&) noexcept = default;
        ~SIMD() = default;
        /* Operators */
        SIMD& operator=(const SIMD&) = default;
        SIMD& operator=(SIMD&&) noexcept = default;
        [[nodiscard]] operator __m128() const noexcept { return toMachine(); }
        [[nodiscard]] operator __m256() const noexcept { return toMachine(); }
        [[nodiscard]] operator __m512() const noexcept { return toMachine(); }
        [[nodiscard]] operator __m128d() const noexcept { return toMachine(); }
        [[nodiscard]] operator __m256d() const noexcept { return toMachine(); }
        [[nodiscard]] operator __m512d() const noexcept { return toMachine(); }
        [[nodiscard]] inline T operator[](int index) const;
        [[nodiscard]] inline SIMD operator+(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator-(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator*(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator*(const T& x) const;
        [[nodiscard]] inline SIMD operator/(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator-() const;
        [[nodiscard]] inline SIMD operator&(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator|(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator^(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator!() const;
        void operator+=(const SIMD& other) { *this = *this + other; }
        void operator-=(const SIMD& other) { *this = *this - other; }
        void operator*=(const SIMD& other) { *this = *this * other; }
        void operator*=(const T& x) { *this = *this * x; }
        void operator/=(const SIMD& other) { *this = *this / other; }
        void operator&=(const SIMD& other) { *this = *this & other; }
        void operator|=(const SIMD& other) { *this = *this | other; }
        void operator^=(const SIMD& other) { *this = *this ^ other; }
        [[nodiscard]] inline auto operator==(const SIMD& other) const;
        [[nodiscard]] inline auto operator>(const SIMD& other) const;
        [[nodiscard]] inline auto operator<(const SIMD& other) const;
        [[nodiscard]] auto operator!=(const SIMD& other) const { return !(*this == other); }
        [[nodiscard]] auto operator>=(const SIMD& other) const { return !(*this < other); }
        [[nodiscard]] auto operator<=(const SIMD& other) const { return !(*this > other); }
        /* Operations */
        inline void load(const T* p);
        inline void load_partial(const T* p, int n);
        inline void store(T* p) const;
        inline void store_partial(T* p, int n) const;
        inline void insert(int index, const T& value);
        template<int... Order> inline SIMD shuffle() const;
        template<int... Order> inline SIMD permute() const;
        template<int... Flags> inline SIMD change_sign() const;

        inline This& cutoff(int count);

        [[nodiscard]] inline T sum() const;
        [[nodiscard]] inline T max() const;
        [[nodiscard]] inline T min() const;
        void swap(SIMD& __restrict other) noexcept { std::swap(*this, other); }
        /* Getters */
        [[nodiscard]] constexpr static size_t size() { return Size; }
        [[nodiscard]] Base& toMachine() noexcept { return *this; }
        [[nodiscard]] const Base& toMachine() const noexcept { return *this; }
        [[nodiscard]] HalfType getLow() const noexcept { return Base::get_low(); }
        [[nodiscard]] HalfType getHigh() const noexcept { return Base::get_high(); }
        [[nodiscard]] This getValue() const noexcept { return *this; }
        [[nodiscard]] auto isPositive() const noexcept { return operator>(This(0)); }
        [[nodiscard]] auto isNegative() const noexcept { return operator<(This(0)); }
        /* Static members */
        template<bool... Flags>
        [[nodiscard]] static SIMD makeSignBits();
        template<RandomGenerator R>
        [[nodiscard]] static SIMD random_uniform();
        [[nodiscard]] static SIMD select(BoolSIMDType flags, const SIMD& x, const SIMD& y);
    private:
        template<int Order, int... Orders>
        [[nodiscard]] constexpr static unsigned int makeShuffleMask(int order);
    };

    template<Scalar T, size_t Size>
    [[nodiscard]] inline SIMD<T, Size> operator*(const T& scalar, const SIMD<T, Size>& packet) {
        return packet * scalar;
    }

    template<Scalar T, size_t Size>
    [[nodiscard]] inline SIMD<T, Size> mul_add(
            const SIMD<T, Size>& a,
            const SIMD<T, Size>& b,
            const SIMD<T, Size>& c);

    template<Scalar T, size_t Size>
    [[nodiscard]] inline SIMD<T, Size> nmul_add(
            const SIMD<T, Size> a,
            const SIMD<T, Size> b,
            const SIMD<T, Size> c);

    template<Scalar T, size_t Size>
    [[nodiscard]] inline SIMD<T, Size> mul_sub(
            const SIMD<T, Size> a,
            const SIMD<T, Size> b,
            const SIMD<T, Size> c);

    template<Scalar T, size_t Size>
    [[nodiscard]] inline SIMD<T, Size> mul_addsub(
            const SIMD<T, Size> a,
            const SIMD<T, Size> b,
            const SIMD<T, Size> c);
}

namespace Physica {
    template<Scalar T, size_t Size>
    class Traits<Core::SIMD<T, Size>> {
        constexpr static bool isFloat32 = T::Option == Core::Float32;
        static_assert(isFloat32 || T::Option == Core::Float64, "[Error]: Unsupported float type");
        static_assert(!T::isComplex, "[Error]: The main template targets on real scalar");
        static_assert(!T::isDifferentiable, "[Error]: The main template targets on plain scalar");
        static_assert(Size % 2 == 0 && Size <= 16, "[Error]: Invalid Size");

        using ValueType = T::ValueType;
        using Size2Type = std::conditional<isFloat32, void, Vec2d>::type;
        using Size4Type = std::conditional<isFloat32, Vec4f, Vec4d>::type;
        using Size8Type = std::conditional<isFloat32, Vec8f, Vec8d>::type;
        using Size16Type = std::conditional<isFloat32, Vec16f, void>::type;
        using Base1 = std::conditional<Size == 2, Size2Type, Size4Type>::type;
        using Base2 = std::conditional<Size == 8, Size8Type, Size16Type>::type;
    public:
        using ScalarType = T;
        using BaseType = std::conditional<Size <= 4, Base1, Base2>::type;
        using BoolSIMDType = Core::BoolSIMD<T, Size>;
    };
}

namespace std {
    #define PacketType Physica::Core::SIMD<T, Size>

    template<Physica::Core::Scalar T, size_t Size>
    inline PacketType max(PacketType a, PacketType b) {
        if constexpr (T::isForwardDiff) {
            using GradPacket = PacketType::GradPacket;
            const auto values = max(a.getValue(), b.getValue());
            return PacketType(values, GradPacket::select(values == a.getValue(), a.getGrad(), b.getGrad()));
        }
        else
            return Physica::max(a.toMachine(), b.toMachine());
    }

    template<Physica::Core::Scalar T, size_t Size>
    inline PacketType min(PacketType a, PacketType b) {
        if constexpr (T::isForwardDiff) {
            using GradPacket = PacketType::GradPacket;
            const auto values = min(a.getValue(), b.getValue());
            return PacketType(values, GradPacket::select(values == a.getValue(), a.getGrad(), b.getGrad()));
        }
        else
            return Physica::min(a.toMachine(), b.toMachine());
    }

    #undef PacketType
}

#include "SIMDImpl/BoolSIMD.h"
#include "SIMDImpl/SIMDImpl.h"
#include "SIMDImpl/BestPacket.h"
#ifdef PHYSICA_CUDA
    #include "SIMDImpl/Half2.h"
#endif
