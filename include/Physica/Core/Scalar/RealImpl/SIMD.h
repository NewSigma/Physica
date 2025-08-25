/*
 * Copyright 2023-2025 Weibo He.
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

#include <format>
#pragma GCC diagnostic push // Necessary when install
#pragma GCC diagnostic ignored "-Wunused-function"
    #include <vectorclass/vectorclass.h>
    #include <vectorclass/vectormath_exp.h>
    #include <vectorclass/vectormath_trig.h>
    #include <vectorclass/vectormath_hyp.h>
#pragma GCC diagnostic pop
#include "Physica/Core/Scalar/Scalar.h"
#include "Physica/Core/Scalar/ScalarImpl/SIMDBase.h"
#include "Physica/PlainStruct.h"

namespace Physica {
    template<Scalar T, int Size> class BoolSIMD;

    template<Scalar T, int Size>
    class SIMD : public SIMDBase<SIMD<T, Size>> {
        constexpr static bool isFloat32 = T::Prec == Float32;
        constexpr static bool isForward = T::isForwardDiff;
        using This = SIMD<T, Size>;
        using Base = SIMDBase<This>;
        using Pack = Traits<This>::Pack;
        using HalfType = std::conditional<sizeof(Pack) * CHAR_BIT != 128, SIMD<T, Size / 2>, PlainStruct<void>>::type;
    public:
        using typename Base::MachineType;
        using Base::isSeparatable;
        using BoolSIMDType = Traits<This>::BoolSIMDType;
    private:
        Pack pack;
    public:
        SIMD() = default;
        explicit SIMD(double x) : pack(x) {}
        explicit SIMD(const Scalar auto& x);
        SIMD(T x, int count);
        SIMD(Scalar auto... args);
        SIMD(MachineType x) : pack(x) {}
        SIMD(Pack value) : pack(value) {}
        SIMD(HalfType a, HalfType b);
        SIMD(const SIMD&) = default;
        SIMD(SIMD&&) noexcept = default;
        ~SIMD() = default;
        /* Operators */
        SIMD& operator=(const SIMD&) = default;
        SIMD& operator=(SIMD&&) noexcept = default;
        [[nodiscard]] operator MachineType() const noexcept { return toMachine(); }
        [[nodiscard]] T operator[](int index) const;
        [[nodiscard]] SIMD operator+(const SIMD& other) const;
        [[nodiscard]] SIMD operator-(const SIMD& other) const;
        [[nodiscard]] SIMD operator*(const SIMD& other) const;
        [[nodiscard]] SIMD operator*(const T& x) const;
        [[nodiscard]] SIMD operator/(const SIMD& other) const;
        [[nodiscard]] SIMD operator-() const;
        [[nodiscard]] SIMD operator&(const SIMD& other) const;
        [[nodiscard]] SIMD operator|(const SIMD& other) const;
        [[nodiscard]] SIMD operator^(const SIMD& other) const;
        [[nodiscard]] SIMD operator!() const;
        void operator+=(const SIMD& other) { *this = *this + other; }
        void operator-=(const SIMD& other) { *this = *this - other; }
        void operator*=(const SIMD& other) { *this = *this * other; }
        void operator*=(const T& x) { *this = *this * x; }
        void operator/=(const SIMD& other) { *this = *this / other; }
        void operator&=(const SIMD& other) { *this = *this & other; }
        void operator|=(const SIMD& other) { *this = *this | other; }
        void operator^=(const SIMD& other) { *this = *this ^ other; }
        [[nodiscard]] auto operator==(const SIMD& other) const;
        [[nodiscard]] auto operator>(const SIMD& other) const;
        [[nodiscard]] auto operator<(const SIMD& other) const;
        [[nodiscard]] auto operator!=(const SIMD& other) const { return !(*this == other); }
        [[nodiscard]] auto operator>=(const SIMD& other) const { return !(*this < other); }
        [[nodiscard]] auto operator<=(const SIMD& other) const { return !(*this > other); }
        /* Operations */
        void load(const T* p);
        void load_partial(const T* p, int n);
        void store(T* p) const;
        void store_partial(T* p, int n) const;
        void insert(int index, const T& value);
        template<int... Order> SIMD shuffle() const;
        template<int... Order> SIMD permute() const;
        template<int... Flags> SIMD change_sign() const;

        This& cutoff(int count);

        [[nodiscard]] T sum() const;
        [[nodiscard]] T max() const;
        [[nodiscard]] T min() const;
        void swap(SIMD& __restrict other) noexcept { std::swap(*this, other); }
        /* Getters */
        [[nodiscard]] Pack& toMachine() noexcept { return pack; }
        [[nodiscard]] const Pack& toMachine() const noexcept { return pack; }
        [[nodiscard]] HalfType getLow() const noexcept { return pack.get_low(); }
        [[nodiscard]] HalfType getHigh() const noexcept { return pack.get_high(); }
        [[nodiscard]] This value() const noexcept { return *this; }
        [[nodiscard]] auto isPositive() const noexcept { return operator>(This(0)); }
        [[nodiscard]] auto isNegative() const noexcept { return operator<(This(0)); }
        [[nodiscard]] BoolSIMDType isFinite() const noexcept;
        /* Static members */
        template<int... Order>
        static SIMD blend(const SIMD& x, const SIMD& y);
        template<bool... Flags>
        [[nodiscard]] static SIMD makeSignBits();
        template<RNG R>
        [[nodiscard]] static SIMD random_uniform();
        [[nodiscard]] static SIMD select(BoolSIMDType flags, const SIMD& x, const SIMD& y);
    private:
        template<int Order, int... Orders>
        [[nodiscard]] constexpr static unsigned int makeShuffleMask(int order);
    };

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<T, Size> operator*(const T& scalar, const SIMD<T, Size>& packet) {
        return packet * scalar;
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<T, Size> mul_add(
            const SIMD<T, Size>& a,
            const SIMD<T, Size>& b,
            const SIMD<T, Size>& c) noexcept;

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<T, Size> nmul_add(
            const SIMD<T, Size> a,
            const SIMD<T, Size> b,
            const SIMD<T, Size> c) noexcept;

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<T, Size> mul_sub(
            const SIMD<T, Size> a,
            const SIMD<T, Size> b,
            const SIMD<T, Size> c) noexcept;

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<T, Size> mul_addsub(
            const SIMD<T, Size> a,
            const SIMD<T, Size> b,
            const SIMD<T, Size> c) noexcept;
}

namespace Physica {
    template<Scalar T, int S>
    class Traits<SIMD<T, S>> {
    public:
        constexpr static int Size = S;
    private:
        constexpr static bool isFloat32 = T::Prec == Float32;
        static_assert(isFloat32 || T::Prec == Float64, "[Error]: Unsupported float type");
        static_assert(!T::isComplex, "[Error]: The main template targets on real scalar");
        static_assert(!T::isDiffable, "[Error]: The main template targets on plain scalar");
        static_assert(Size % 2 == 0 && Size <= 16, "[Error]: Invalid Size");

        using Size2Type = std::conditional<isFloat32, void, Vec2d>::type;
        using Size4Type = std::conditional<isFloat32, Vec4f, Vec4d>::type;
        using Size8Type = std::conditional<isFloat32, Vec8f, Vec8d>::type;
        using Size16Type = std::conditional<isFloat32, Vec16f, void>::type;
        using Pack1 = std::conditional<Size == 2, Size2Type, Size4Type>::type;
        using Pack2 = std::conditional<Size == 8, Size8Type, Size16Type>::type;

        template<class, int I = 0> struct MachineTypeHelper; // We have to use specialization to avoid discarding type attribute
        template<int I> struct MachineTypeHelper<Vec4f, I> { using Type = __m128; };
        template<int I> struct MachineTypeHelper<Vec8f, I> { using Type = __m256; };
        template<int I> struct MachineTypeHelper<Vec16f, I> { using Type = __m512; };
        template<int I> struct MachineTypeHelper<Vec2d, I> { using Type = __m128d; };
        template<int I> struct MachineTypeHelper<Vec4d, I> { using Type = __m256d; };
        template<int I> struct MachineTypeHelper<Vec8d, I> { using Type = __m512d; };
    public:
        using Pack = std::conditional<Size <= 4, Pack1, Pack2>::type;

        using ScalarType = T;
        using ValueType = SIMD<T, Size>;
        using GradType = void;
        using RealType = ValueType;
        using FullRealType = ValueType;
        using MachineType = MachineTypeHelper<Pack>::Type;
        using BoolSIMDType = BoolSIMD<T, Size>;

        constexpr static bool isSeparatable = sizeof(MachineType) * CHAR_BIT != 128;
        static_assert(!std::is_same<Pack, void>::value, "[Error]: Bad packet");
    };
}

namespace std {
#define PacketType Physica::SIMD<T, Size>

    template<Physica::Scalar T, int Size>
    PacketType max(PacketType a, PacketType b) {
        static_assert(!T::isComplex, "[Error]: Compare between complex number is ill defined");
        if constexpr (T::isForwardDiff) {
            using GradPacket = PacketType::GradType;
            const auto values = max(a.value(), b.value());
            return PacketType(values, GradPacket::select(values == a.value(), a.grad(), b.grad()));
        }
        else
            return Physica::max(a.toMachine(), b.toMachine());
    }

    template<Physica::Scalar T, int Size>
    PacketType min(PacketType a, PacketType b) {
        static_assert(!T::isComplex, "[Error]: Compare between complex number is ill defined");
        if constexpr (T::isForwardDiff) {
            using GradPacket = PacketType::GradType;
            const auto values = min(a.value(), b.value());
            return PacketType(values, GradPacket::select(values == a.value(), a.grad(), b.grad()));
        }
        else
            return Physica::min(a.toMachine(), b.toMachine());
    }

    template<Physica::Scalar T, int Size>
    struct formatter<PacketType, char> {
        constexpr auto parse(std::format_parse_context& ctx) {
            return ctx.begin();
        }

        auto format(const PacketType& obj, std::format_context& ctx) const {
            for (size_t i = 0; i < Size - 1; ++i)
                std::format_to(ctx.out(), "{} ", obj[i]);
            return std::format_to(ctx.out(), "{}", obj[Size - 1]);
        }
    };

#undef PacketType
}

#include "SIMDImpl/BestPacket.h"
#include "SIMDImpl/BoolSIMD.h"
#include "SIMDImpl/SIMDImpl.h"
#ifdef PHYSICA_CUDA
    #include "SIMDImpl/Half2.h"
#endif
