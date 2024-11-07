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
    template<class ScalarType, int Order> class DiffTracer;

    template<class T, int Order, size_t Length>
    class BestPacket<Diff<T, DiffMode::Forward, Order>, Length> {
        using ScalarType = Diff<T, DiffMode::Forward, Order>;
    public:
        constexpr static size_t Size = BestPacket<T, Length>::Size;
        using Type = typename std::conditional<Size == 1, ScalarType, SIMD<ScalarType, Size>>::type;
    };

    template<class T, int Order, size_t Length>
    class device_obj<BestPacket<Diff<T, DiffMode::Forward, Order>, Length>> {
    public:
        constexpr static size_t Size = 1;
        using Type = Diff<T, DiffMode::Forward, Order>;
    };

    template<class T, int Order, size_t Size>
    class SIMD<Diff<T, DiffMode::Forward, Order>, Size> {
        using ScalarType = Diff<T, DiffMode::Forward, Order>;
        using This = SIMD<ScalarType, Size>;
        using PtrTy = typename ScalarType::PtrTy;
        using ConstPtrTy = typename ScalarType::ConstPtrTy;
        using GradType = typename ScalarType::GradType;
    public:
        using ValuePacket = SIMD<T, Size>;
        using GradPacket = SIMD<GradType, Size>;
        using BoolSIMDType = typename ValuePacket::BoolSIMDType;
        using HalfType = typename std::conditional<sizeof(ValuePacket) * CHAR_BIT != 128, SIMD<ScalarType, Size / 2>, PlainStruct<void>>::type;

        ValuePacket values;
        GradPacket grads;
    public:
        SIMD() = default;
        explicit SIMD(int x);
        explicit SIMD(double x);
        explicit SIMD(ScalarType x);
        SIMD(ScalarType x, int count);
        explicit SIMD(ValuePacket values_);
        SIMD(ValuePacket values_, GradPacket grads_);
        template<int OtherOrder>
        explicit SIMD(const SIMD<Diff<T, DiffMode::Forward, OtherOrder>, Size>& other);
        SIMD(const SIMD&) = default;
        SIMD(SIMD&&) noexcept = default;
        ~SIMD() = default;
        /* Operators */
        SIMD& operator=(const SIMD&) = default;
        SIMD& operator=(SIMD&&) noexcept = default;
        [[nodiscard]] explicit operator ValuePacket() const noexcept { return values; }
        [[nodiscard]] inline ScalarType operator[](int index) const;
        [[nodiscard]] inline SIMD operator+(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator-(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator*(const SIMD& x) const;
        [[nodiscard]] inline SIMD operator*(const ScalarType& x) const;
        //[[nodiscard]] inline SIMD operator*(const T& x) const;
        [[nodiscard]] inline SIMD operator/(const SIMD& x) const;
        [[nodiscard]] inline SIMD operator-() const;
        void operator+=(const SIMD& x) { *this = *this + x; }
        void operator-=(const SIMD& x) { *this = *this - x; }
        void operator*=(const SIMD& x) { *this = *this * x; }
        void operator/=(const SIMD& x) { *this = *this / x; }
        /* Operations */
        inline void load(ConstPtrTy p);
        inline void load_partial(ConstPtrTy p, int n);
        inline void store(PtrTy p) const;
        inline void store_partial(PtrTy p, int n) const;

        inline This& cutoff(int count);

        [[nodiscard]] inline ScalarType sum() const;
        [[nodiscard]] inline ScalarType max() const;
        [[nodiscard]] inline ScalarType min() const;
        void swap(SIMD& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] constexpr static size_t size() { return Size; }
        [[nodiscard]] ValuePacket getValue() const noexcept { return values; }
        [[nodiscard]] GradPacket getGrad() const noexcept { return grads; }
        [[nodiscard]] HalfType getLow() const noexcept { return HalfType(values.getLow(), grads.getLow()); }
        [[nodiscard]] HalfType getHigh() const noexcept { return HalfType(values.getHigh(), grads.getHigh()); }
        /* Static members */
        [[nodiscard]] static SIMD select(BoolSIMDType flags, const SIMD& x, const SIMD& y);
    };

    template<class PlainScalar, size_t Size>
    class SIMD<Diff<PlainScalar, DiffMode::Reverse, 1>, Size> : public SIMD<PlainScalar, Size> {
        static_assert(!PlainScalar::isDifferentiable, "[Error]: Invalid template param");
        using ScalarType = Diff<PlainScalar, DiffMode::Reverse, 1>;
        using This = SIMD<ScalarType, Size>;
        using Base = SIMD<PlainScalar, Size>;
        using BoolSIMDType = typename Traits<This>::BoolSIMDType;
    public:
        using PlainPacket = Base;
        using TracerType = DiffTracer<PlainScalar, 1>;
    private:
        ScalarType headNode;
    public:
        SIMD() = default;
        explicit SIMD(int i) : SIMD(PlainScalar(i)) {}
        explicit SIMD(PlainScalar s);
        explicit SIMD(ScalarType s);
        SIMD(Base base, ScalarType headNode_);
        SIMD(const SIMD&) = default;
        SIMD(SIMD&&) noexcept = default;
        ~SIMD() = default;
        /* Operators */
        SIMD& operator=(const SIMD&) = default;
        SIMD& operator=(SIMD&&) noexcept = default;
        [[nodiscard]] inline ScalarType operator[](int i) const;
        [[nodiscard]] inline SIMD operator+(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator-(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator*(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator*(const ScalarType& v) const;
        [[nodiscard]] inline SIMD operator/(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator-() const;
        void operator+=(const SIMD& x) { *this = *this + x; }
        void operator-=(const SIMD& x) { *this = *this - x; }
        void operator*=(const SIMD& x) { *this = *this * x; }
        void operator/=(const SIMD& x) { *this = *this / x; }
        using Base::operator>;
        using Base::operator<;
        using Base::operator>=;
        using Base::operator<=;
        /* Operations */
        inline void load(const ScalarType* p);
        inline void load_partial(const ScalarType* p, int n);
        inline void store(ScalarType* p) const;
        inline void store_partial(ScalarType* p, int n) const;
        [[nodiscard]] inline ScalarType sum() const;
        inline void swap(SIMD& __restrict obj) noexcept;
        /* Getters */
        using Base::size;
        using Base::getImpl;
        [[nodiscard]] ScalarType getHeadNode() const noexcept { return headNode; }
        [[nodiscard]] PlainScalar* value_ptr() const noexcept { return headNode.value_ptr(); }
        [[nodiscard]] PlainScalar* grad_ptr() const noexcept { return headNode.grad_ptr(); }
    private:
        using Base::insert; //Insert a scalar may lead to incontineous memory, which harms performance
        [[nodiscard]] static bool checkContinuous(const ScalarType* p, int n);
    };

    template<class T, DiffMode Mode, int Order, size_t Size>
    [[nodiscard]] inline SIMD<Diff<T, Mode, Order>, Size> mul_add(
            const SIMD<Diff<T, Mode, Order>, Size>& a,
            const SIMD<Diff<T, Mode, Order>, Size>& b,
            const SIMD<Diff<T, Mode, Order>, Size>& c);
}

namespace Physica {
    template<class T, Core::DiffMode Mode, int Order, size_t Size>
    class Traits<Core::SIMD<Core::Diff<T, Mode, Order>, Size>> : public Traits<Core::SIMD<T, Size>> {
    public:
        using ScalarType = Core::Diff<T, Mode, Order>;
    };
}

#include "SIMDImpl/SIMDImpl.h"
#include "SIMDImpl/Math.h"
