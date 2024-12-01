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

namespace Physica::Core {
    template<Scalar T, size_t Length>
    class BestPacket<Complex<T>, Length> {
        using RealPacket = BestPacket<T, Length * 2>;
        using ScalarType = Complex<T>;
    public:
        constexpr static size_t Size = RealPacket::Size / 2;
        using Type = std::conditional<Size == 1, ScalarType, SIMD<ScalarType, Size>>::type;
    };

    template<Scalar T, size_t Length>
    class device_obj<BestPacket<Complex<T>, Length>> {
    public:
        constexpr static size_t Size = 1;
        using Type = Complex<T>;
    };

    template<Scalar T, size_t Size>
    class SIMD<Complex<T>, Size> : public SIMDBase<SIMD<Complex<T>, Size>>, private SIMD<T, Size * 2> {
        using ScalarType = Complex<T>;
        using This = SIMD<ScalarType, Size>;
        using Base = SIMDBase<This>;
        using RealBase = SIMD<T, Size * 2>;
        using MachineType = ScalarType::MachineType;
    public:
        using typename Base::RealType;
        using typename Base::FullRealType;
        using Base::isSeparatable;
        using BoolSIMDType = FullRealType::BoolSIMDType;
        using PlainPacket = This;
    private:
        using HalfType = std::conditional<sizeof(FullRealType) * CHAR_BIT != 128, SIMD<Complex<T>, Size / 2>, PlainStruct<void>>::type;

        using FullRealPair = std::pair<FullRealType, FullRealType>;
    public:
        SIMD() = default;
        explicit SIMD(int x);
        explicit SIMD(ScalarType x);
        SIMD(const SIMD&) = default;
        SIMD(SIMD&&) noexcept = default;
        ~SIMD() = default;
        /* Operators */
        SIMD& operator=(const SIMD&) = default;
        SIMD& operator=(SIMD&&) noexcept = default;
        [[nodiscard]] inline ScalarType operator[](int index) const;
        [[nodiscard]] inline SIMD operator+(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator-(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator*(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator*(const ScalarType& x) const;
        [[nodiscard]] inline SIMD operator*(const T& x) const;
        [[nodiscard]] inline SIMD operator/(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator-() const;
        void operator+=(const SIMD& other) { *this = *this + other; }
        void operator-=(const SIMD& other) { *this = *this - other; }
        void operator*=(const SIMD& other) { *this = *this * other; }
        void operator/=(const SIMD& other) { *this = *this / other; }
        /* Operations */
        inline void load(const ScalarType* p);
        inline void load_partial(const ScalarType* p, int n);
        inline void store(ScalarType* p) const;
        inline void store_partial(ScalarType* p, int n) const;

        inline FullRealPair makeFullRealImag() const noexcept;
        using Base::swapRealImag;
        using Base::permRealImag;
        using Base::squaredNorm;

        //inline void insert(int index, const ScalarType& value);
        [[nodiscard]] inline ScalarType sum() const;
        void swap(SIMD& __restrict other) noexcept { std::swap(*this, other); }
        /* Getters */
        using Base::getValue;
        [[nodiscard]] constexpr static size_t size() { return Size; }
        [[nodiscard]] FullRealType asReal() const noexcept { return FullRealType::toMachine(); }
        [[nodiscard]] HalfType getLow() const noexcept { return HalfType::asComplex(FullRealType::getLow()); }
        [[nodiscard]] HalfType getHigh() const noexcept { return HalfType::asComplex(FullRealType::getHigh()); }
        [[nodiscard]] inline RealType real() const noexcept;
        [[nodiscard]] inline RealType imag() const noexcept;
        /* Static members */
        template<RandomGenerator R>
        [[nodiscard]] static SIMD random_uniform() { return asComplex(FullRealType::template random_uniform<R>()); }
        [[nodiscard]] static SIMD asComplex(FullRealType reals);
    };
}

namespace Physica {
    template<Scalar T, size_t S>
    class Traits<SIMD<Complex<T>, S>> {
    public:
        constexpr static int Size = S;

        using ScalarType = Complex<T>;
        using ValueType = SIMD<Complex<T>, S>;
        using FullRealType = SIMD<T, Size * 2>;
        using RealType = std::conditional<FullRealType::isSeparatable, SIMD<T, Size>, FullRealType>::type;
        using MachineType = FullRealType::MachineType;
        constexpr static bool isSeparatable = FullRealType::isSeparatable;
    };
}

#include "SIMDImpl/SIMDImpl.h"
#include "SIMDImpl/Math.h"
