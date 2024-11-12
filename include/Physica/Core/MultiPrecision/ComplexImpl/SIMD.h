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
    template<class T, size_t Length>
    class BestPacket<Complex<T>, Length> {
        using RealPacket = BestPacket<T, Length * 2>;
        using ScalarType = Complex<T>;
    public:
        constexpr static size_t Size = RealPacket::Size / 2;
        using Type = typename std::conditional<Size == 1, ScalarType, SIMD<ScalarType, Size>>::type;
    };

    template<class ScalarType, size_t Length>
    class device_obj<BestPacket<Complex<ScalarType>, Length>> {
    public:
        constexpr static size_t Size = 1;
        using Type = Complex<ScalarType>;
    };

    template<class T, size_t Size>
    class SIMD<Complex<T>, Size> : private SIMD<T, Size * 2> {
        using ScalarType = Complex<T>;
        using This = SIMD<ScalarType, Size>;
        using Base = SIMD<T, Size * 2>;
        using RealType = SIMD<T, Size>;
        using AsRealRtnTy = Base;
        using MachineType = typename ScalarType::MachineType;
        using HalfType = typename std::conditional<sizeof(Base) * CHAR_BIT != 128, SIMD<Complex<T>, Size / 2>, PlainStruct<void>>::type;
        using Base::isSeparatable;
    public:
        using BoolSIMDType = typename Base::BoolSIMDType;
        using PlainPacket = This;
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
        //[[nodiscard]] inline SIMD operator/(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator-() const;
        void operator+=(const SIMD& other) { *this = *this + other; }
        void operator-=(const SIMD& other) { *this = *this - other; }
        void operator*=(const SIMD& other) { *this = *this * other; }
        //void operator/=(const SIMD& other) { *this = *this / other; }
        /* Operations */
        inline void load(const ScalarType* p);
        inline void load_partial(const ScalarType* p, int n);
        inline void store(ScalarType* p) const;
        inline void store_partial(ScalarType* p, int n) const;

        [[nodiscard]] inline Base swapRealImag() const noexcept;
        [[nodiscard]] inline Base permRealImag() const noexcept;

        //inline void insert(int index, const ScalarType& value);
        [[nodiscard]] inline ScalarType sum() const;
        void swap(SIMD& __restrict other) noexcept { std::swap(*this, other); }
        /* Getters */
        [[nodiscard]] constexpr static size_t size() { return Size; }
        [[nodiscard]] AsRealRtnTy asReal() const noexcept { return Base::toMachine(); }
        [[nodiscard]] HalfType getLow() const noexcept { return HalfType::asComplex(Base::getLow()); }
        [[nodiscard]] HalfType getHigh() const noexcept { return HalfType::asComplex(Base::getHigh()); }
        [[nodiscard]] inline RealType real() const noexcept;
        [[nodiscard]] inline RealType imag() const noexcept;
        /* Static members */
        template<class RandomType>
        [[nodiscard]] static SIMD random_uniform(RandomType& gen) { return asComplex(Base::random_uniform(gen)); }
        [[nodiscard]] static SIMD asComplex(AsRealRtnTy reals);
    };
}

namespace Physica {
    template<class T, size_t Size>
    class Traits<Core::SIMD<Core::Complex<T>, Size>> {
    public:
        using ScalarType = Core::Complex<T>;
    };
}

#include "SIMDImpl/SIMDImpl.h"
#include "SIMDImpl/Math.h"
