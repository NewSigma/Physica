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

#include "../Diff.h"

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order, size_t Length>
    class BestPacket<Diff<T, Mode, Order>, Length> {
        using ScalarType = Diff<T, Mode, Order>;
    public:
        constexpr static size_t Size = BestPacket<T, Length>::Size;
        using Type = std::conditional<Size <= 1, ScalarType, SIMD<ScalarType, Size>>::type;
    };

    template<Scalar T, DiffMode Mode, int Order, size_t Length>
    class device_obj<BestPacket<Diff<T, Mode, Order>, Length>> {
    public:
        constexpr static size_t Size = 1;
        using Type = Diff<T, Mode, Order>;
    };

    template<DiffMode Mode, int Order, size_t Length>
    class device_obj<BestPacket<Diff<float16, Mode, Order>, Length>> {
        using ScalarType = std::conditional<Mode == DiffMode::Forward, Diff<float16, Mode, Order>, float16>::type;
    public:
        constexpr static size_t Size = Length == 1 ? 1 : 2;
        using Type = SIMD<ScalarType, Size>;
    };

    template<Scalar T, DiffMode Mode, int Order, size_t Size>
    class SIMD<Diff<T, Mode, Order>, Size> : public SIMDBase<SIMD<Diff<T, Mode, Order>, Size>> {
        using This = SIMD<Diff<T, Mode, Order>, Size>;
        using Base = SIMDBase<This>;
        using RealType = SIMD<Diff<typename T::RealType, Mode, Order>, Size>;
        using FullRealType = SIMD<Diff<typename T::RealType, Mode, Order>, Size * (T::isComplex ? 2 : 1)>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using typename Base::GradType;
        using typename Base::BoolSIMDType;
        using PtrTy = ScalarType::PtrTy;
        using ConstPtrTy = ScalarType::ConstPtrTy;
        using HalfType = std::conditional<sizeof(ValueType) * CHAR_BIT != 128, SIMD<ScalarType, Size / 2>, PlainStruct<void>>::type;
    private:
        ValueType values;
        GradType grads;
    public:
        SIMD() = default;
        explicit SIMD(int x);
        explicit SIMD(double x);
        explicit SIMD(Scalar auto x);
        SIMD(ScalarType x, int count);
        explicit SIMD(ValueType values_);
        SIMD(ValueType values_, GradType grads_);
        SIMD(HalfType a, HalfType b);
        template<int OtherOrder>
        explicit SIMD(const SIMD<Diff<T, Mode, OtherOrder>, Size>& other);
        SIMD(const SIMD&) = default;
        SIMD(SIMD&&) noexcept = default;
        ~SIMD() = default;
        /* Operators */
        SIMD& operator=(const SIMD&) = default;
        SIMD& operator=(SIMD&&) noexcept = default;
        [[nodiscard]] explicit operator ValueType() const noexcept { return values; }
        [[nodiscard]] inline ScalarType operator[](int index) const;
        [[nodiscard]] inline SIMD operator+(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator-(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator*(const SIMD& x) const;
        [[nodiscard]] inline SIMD operator*(const ScalarType& x) const;
        [[nodiscard]] inline SIMD operator*(const Scalar auto& x) const;
        [[nodiscard]] inline SIMD operator/(const SIMD& x) const;
        [[nodiscard]] inline SIMD operator-() const;
        void operator+=(const SIMD& x) { *this = *this + x; }
        void operator-=(const SIMD& x) { *this = *this - x; }
        void operator*=(const SIMD& x) { *this = *this * x; }
        void operator/=(const SIMD& x) { *this = *this / x; }
        /* Operations */
        ValueType reverse(GradType grad = 1) const noexcept;

        inline void load(ConstPtrTy p);
        inline void load_partial(ConstPtrTy p, int n);
        inline void store(PtrTy p) const;
        inline void store_partial(PtrTy p, int n) const;

        inline This& cutoff(int count);
        [[nodiscard]] inline FullRealType swapRealImag() const noexcept;
        [[nodiscard]] inline FullRealType permRealImag() const noexcept;

        [[nodiscard]] inline auto sum() const;
        [[nodiscard]] inline auto max() const;
        [[nodiscard]] inline auto min() const;
        void swap(SIMD& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] FullRealType asReal() const noexcept;
        [[nodiscard]] ValueType& value() noexcept { return values; }
        [[nodiscard]] const ValueType& value() const noexcept { return values; }
        [[nodiscard]] GradType& grad() noexcept { return grads; }
        [[nodiscard]] const GradType& grad() const noexcept { return grads; }
        [[nodiscard]] HalfType getLow() const noexcept { return HalfType(values.getLow(), grads.getLow()); }
        [[nodiscard]] HalfType getHigh() const noexcept { return HalfType(values.getHigh(), grads.getHigh()); }
        [[nodiscard]] inline RealType real() const noexcept;
        [[nodiscard]] inline RealType imag() const noexcept;
        /* Static members */
        [[nodiscard]] static SIMD select(BoolSIMDType flags, const SIMD& x, const SIMD& y);
        [[nodiscard]] static SIMD asComplex(const FullRealType& reals);
    };

    template<Scalar T, DiffMode Mode, int Order, size_t Size>
    [[nodiscard]] inline auto mul_add(
            const SIMD<Diff<T, Mode, Order>, Size>& a,
            const SIMD<Diff<T, Mode, Order>, Size>& b,
            const SIMD<Diff<T, Mode, Order>, Size>& c);
}

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order, size_t Size>
    class Traits<SIMD<Diff<T, Mode, Order>, Size>> : public Traits<SIMD<T, Size>> {
    public:
        using ScalarType = Diff<T, Mode, Order>;
        using GradType = SIMD<typename ScalarType::GradType, Size>;
    };
}

#include "SIMDImpl/SIMDImpl.h"
#include "SIMDImpl/Math.h"
