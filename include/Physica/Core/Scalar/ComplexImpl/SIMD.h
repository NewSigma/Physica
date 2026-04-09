/*
 * Copyright 2024-2026 Weibo He.
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

#include "Physica/Core/Scalar/Complex.h"

namespace Physica {
    template<Scalar T, size_t Length>
    class BestPacket<Complex<T>, Length> {
        using RealPacket = BestPacket<T, Length * 2>;
        using ScalarType = Complex<T>;
    public:
        constexpr static int Size = RealPacket::Size / 2;
        using Type = std::conditional<Size <= 1, ScalarType, SIMD<ScalarType, Size>>::type;
    };

    template<Scalar T, size_t Length>
    class device_obj<BestPacket<Complex<T>, Length>> {
    public:
        constexpr static int Size = 1;
        using Type = Complex<T>;
    };

    template<Scalar T, int Size>
    class SIMD<Complex<T>, Size> : public SIMDBase<SIMD<Complex<T>, Size>> {
        using This = SIMD<Complex<T>, Size>;
        using Base = SIMDBase<This>;
        using MachineType = Complex<T>::MachineType;
    public:
        using typename Base::ScalarType;
        using typename Base::RealType;
        using typename Base::FullRealType;
        using typename Base::BoolSIMDType;
        using Base::isSeparatable;
    private:
        using HalfType = std::conditional<sizeof(FullRealType) * CHAR_BIT != 128, SIMD<Complex<T>, Size / 2>, PlainStruct<void>>::type;
        using FullRealPair = std::pair<FullRealType, FullRealType>;

        FullRealType storage;
    public:
        SIMD() = default;
        explicit SIMD(int x);
        explicit SIMD(double x);
        explicit SIMD(ScalarType x);
        SIMD(ScalarType x, int count);
        SIMD(const SIMD&) = default;
        SIMD(SIMD&&) noexcept = default;
        ~SIMD() = default;
        /* Operators */
        SIMD& operator=(const SIMD&) = default;
        SIMD& operator=(SIMD&&) noexcept = default;
        [[nodiscard]] ScalarType operator[](int index) const;
        [[nodiscard]] auto operator+(Packet auto x) const;
        [[nodiscard]] auto operator-(Packet auto x) const;
        [[nodiscard]] auto operator*(Packet auto x) const;
        [[nodiscard]] auto operator*(const Scalar auto& x) const;
        [[nodiscard]] auto operator/(Packet auto x) const;
        [[nodiscard]] SIMD operator-() const;
        void operator+=(SIMD other) { *this = *this + other; }
        void operator-=(SIMD other) { *this = *this - other; }
        void operator*=(SIMD other) { *this = *this * other; }
        void operator/=(SIMD other) { *this = *this / other; }
        /* Operations */
        void load(const ScalarType* p) noexcept;
        void load(const ScalarType* p, int n) noexcept;
        void store(ScalarType* p) const noexcept;
        void store(ScalarType* p, int n) const noexcept;

        This& cutoff(int count);

        [[nodiscard]] This conjugate() const noexcept;
        [[nodiscard]] FullRealPair makeFullRealImag() const noexcept;
        using Base::swapRealImag;
        using Base::gatherRealImag;
        using Base::squaredNorm;

        //void insert(int index, const ScalarType& value);
        [[nodiscard]] ScalarType sum() const;
        void swap(SIMD& __restrict other) noexcept { storage.swap(other.storage); }
        /* Getters */
        using Base::size;
        using Base::value;
        [[nodiscard]] auto isZero() const noexcept;
        [[nodiscard]] auto isFinite() const noexcept;
        [[nodiscard]] auto isSubNormal() const noexcept;
        [[nodiscard]] FullRealType asReal() const noexcept { return storage; }
        [[nodiscard]] HalfType getLow() const noexcept { return HalfType::asComplex(storage.getLow()); }
        [[nodiscard]] HalfType getHigh() const noexcept { return HalfType::asComplex(storage.getHigh()); }
        [[nodiscard]] RealType real() const noexcept;
        [[nodiscard]] RealType imag() const noexcept;
        /* Static members */
        [[nodiscard]] static SIMD zeros() noexcept;
        template<RNG R>
        [[nodiscard]] static SIMD random_uniform() { return asComplex(FullRealType::template random_uniform<R>()); }
        [[nodiscard, gnu::always_inline]] static SIMD asComplex(FullRealType storage) noexcept;
    };
}

namespace Physica {
    template<Scalar T, int S>
    class Traits<SIMD<Complex<T>, S>> {
    public:
        constexpr static int Size = S;

        using ScalarType = Complex<T>;
        using ValueType = SIMD<Complex<T>, S>;
        using GradType = void;
        using FullRealType = SIMD<T, Size * 2>;
        using RealType = std::conditional<FullRealType::isSeparatable, SIMD<T, Size>, FullRealType>::type;
        using BoolSIMDType = FullRealType::BoolSIMDType;
        using MachineType = FullRealType::MachineType;
        constexpr static bool isSeparatable = FullRealType::isSeparatable;
    };
}

#include "SIMDImpl/SIMDImpl.h"
#include "SIMDImpl/Math.h"
