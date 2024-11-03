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
    public:
        constexpr static size_t Size = 1;
        using Type = Diff<T, DiffMode::Forward, Order>;
    };

    template<class T, int Order, size_t Length>
    class device_obj<BestPacket<Diff<T, DiffMode::Forward, Order>, Length>> {
    public:
        constexpr static size_t Size = 1;
        using Type = Diff<T, DiffMode::Forward, Order>;
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
        void operator+=(const SIMD& other) { *this = *this + other; }
        void operator-=(const SIMD& other) { *this = *this - other; }
        void operator*=(const SIMD& other) { *this = *this * other; }
        void operator/=(const SIMD& other) { *this = *this / other; }
        using Base::operator>;
        using Base::operator<;
        using Base::operator>=;
        using Base::operator<=;
        /* Operations */
        inline void load(const ScalarType* p);
        inline void load_partial(int n, const ScalarType* p);
        inline void store(ScalarType* p) const;
        inline void store_partial(int n, ScalarType* p) const;
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
        [[nodiscard]] static bool checkContinuous(int n, const ScalarType* p);
    };
}

#include "SIMDImpl.h"
