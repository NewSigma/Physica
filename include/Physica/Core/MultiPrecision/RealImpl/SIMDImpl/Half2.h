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

#include <Physica/Core/Utils/Unreachable.h>

namespace Physica::Core {
    template<class ScalarType, size_t Length> class BestPacket;

    template<size_t Length>
    class BestPacket<float16, Length> {
        using ScalarType = float16;
    public:
        constexpr static size_t Size = Length == 1 ? 1 : 2;
        using Type = typename std::conditional<Length == 1, float16, SIMD<float16, 2>>::type;
    };

    template<size_t Length>
    class device_obj<BestPacket<float16, Length>> : public BestPacket<float16, Length> {};

    template<>
    class SIMD<Real<Float16>, 2> : private __half2 {
        constexpr static size_t Size = 2;
        using ScalarType = Real<Float16>;
        using This = SIMD<ScalarType, Size>;
        using Base = __half2;
    public:
        using Base::Base;
        SIMD() = default;
        __host__ __device__ inline explicit SIMD(ScalarType v);
        __host__ __device__ inline SIMD(ScalarType l, ScalarType h);
        __host__ __device__ inline explicit SIMD(Base value);
        SIMD(const SIMD&) = default;
        SIMD(SIMD&&) noexcept = default;
        ~SIMD() = default;
        /* Operators */
        SIMD& operator=(const SIMD&) = default;
        SIMD& operator=(SIMD&&) noexcept = default;
        [[nodiscard]] __host__ __device__ inline ScalarType operator[](int index) const noexcept;
        [[nodiscard]] __host__ __device__ inline SIMD operator+(const SIMD& other) const noexcept;
        [[nodiscard]] __host__ __device__ inline SIMD operator-(const SIMD& other) const noexcept;
        [[nodiscard]] __host__ __device__ inline SIMD operator*(const SIMD& other) const noexcept;
        [[nodiscard]] __host__ __device__ inline SIMD operator/(const SIMD& other) const noexcept;
        [[nodiscard]] __host__ __device__ inline SIMD operator-() const noexcept;
        __host__ __device__ void operator+=(const SIMD& other) noexcept { *this = *this + other; }
        __host__ __device__ void operator-=(const SIMD& other) noexcept { *this = *this - other; }
        __host__ __device__ void operator*=(const SIMD& other) noexcept { *this = *this * other; }
        __host__ __device__ void operator/=(const SIMD& other) noexcept { *this = *this / other; }
        //[[nodiscard]] inline BoolSIMDType operator>(const SIMD& other) const;
        //[[nodiscard]] inline BoolSIMDType operator<(const SIMD& other) const;
        //[[nodiscard]] inline BoolSIMDType operator>=(const SIMD& other) const { return !(*this < other); }
        //[[nodiscard]] inline BoolSIMDType operator<=(const SIMD& other) const { return !(*this > other); }
        /* Operations */
        __host__ __device__ inline void load(const ScalarType* p);
        __host__ __device__ inline void load_partial(const ScalarType* p, int n);
        __host__ __device__ inline void store(ScalarType* p) const;
        __host__ __device__ inline void store_partial(ScalarType* p, int n) const;
        //inline void insert(int index, const ScalarType& value);
        [[nodiscard]] __host__ __device__ inline ScalarType sum() const noexcept;
        //[[nodiscard]] inline ScalarType max() const;
        //[[nodiscard]] inline ScalarType min() const;
        __host__ __device__ void swap(SIMD& __restrict other) noexcept { std::swap(*this, other); }
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr static size_t size() { return Size; }
        [[nodiscard]] __host__ __device__ Base& toMachine() noexcept { return *this; }
        [[nodiscard]] __host__ __device__ const Base& toMachine() const noexcept { return *this; }
    };

    __host__ __device__ inline SIMD<Real<Float16>, 2>::SIMD(ScalarType v) : Base(__half2half2(v.toMachine())) {}

    __host__ __device__ inline SIMD<Real<Float16>, 2>::SIMD(ScalarType l, ScalarType h)
            : Base(make_half2(l.toMachine(), h.toMachine())) {}

    __host__ __device__ inline SIMD<Real<Float16>, 2>::SIMD(Base value) : Base(value) {}

    __host__ __device__ inline Real<Float16> SIMD<Real<Float16>, 2>::operator[](int index) const noexcept {
        if (index == 0)
            return __low2half(*this);
        return __high2half(*this);
    }

    __host__ __device__ inline SIMD<Real<Float16>, 2> SIMD<Real<Float16>, 2>::operator+(const SIMD& other) const noexcept {
        return This(__hadd2(toMachine(), other.toMachine()));
    }

    __host__ __device__ inline SIMD<Real<Float16>, 2> SIMD<Real<Float16>, 2>::operator-(const SIMD& other) const noexcept {
        return This(__hsub2(toMachine(), other.toMachine()));
    }

    __host__ __device__ inline SIMD<Real<Float16>, 2> SIMD<Real<Float16>, 2>::operator*(const SIMD& other) const noexcept {
        return This(__hmul2(toMachine(), other.toMachine()));
    }

    __host__ __device__ inline SIMD<Real<Float16>, 2> SIMD<Real<Float16>, 2>::operator/(const SIMD& other) const noexcept {
        return This(__h2div(toMachine(), other.toMachine()));
    }

    __host__ __device__ inline SIMD<Real<Float16>, 2> SIMD<Real<Float16>, 2>::operator-() const noexcept {
        return This(__hneg2(toMachine()));
    }

    __host__ __device__ inline void SIMD<Real<Float16>, 2>::load(const ScalarType* p) {
        *reinterpret_cast<uint32_t*>(this) = *reinterpret_cast<const uint32_t*>(p);
    }

    __host__ __device__ inline void SIMD<Real<Float16>, 2>::load_partial(const ScalarType* p, int n) {
        if (n == 1)
            (*this) = SIMD(*p, 0);
        else if (n == 2)
            load(p);
    }

    __host__ __device__ inline void SIMD<Real<Float16>, 2>::store(ScalarType* p) const {
        *reinterpret_cast<uint32_t*>(p) = *reinterpret_cast<const uint32_t*>(this);
    }

    __host__ __device__ inline void SIMD<Real<Float16>, 2>::store_partial(ScalarType* p, int n) const {
        if (n == 1)
            *p = operator[](0);
        else if (n == 2)
            store(p);
    }

    __host__ __device__ inline Real<Float16> SIMD<Real<Float16>, 2>::sum() const noexcept {
        return operator[](0) + operator[](1);
    }

    [[nodiscard]] __host__ __device__ inline SIMD<Real<Float16>, 2> mul_add(
            const SIMD<Real<Float16>, 2>& a, const SIMD<Real<Float16>, 2>& b, const SIMD<Real<Float16>, 2>& c) noexcept {
    #ifdef __CUDA_ARCH__
        return SIMD<Real<Float16>, 2>(__hfma2(a.toMachine(), b.toMachine(), c.toMachine()));
    #else
        return a * b + c;
    #endif
    }
}

namespace Physica {
    #define T Core::Real<Core::Float16>

    template<>
    class Traits<Core::SIMD<T, 2>> {
    public:
        using ScalarType = T;
        using BaseType = __half2;
    };

    #undef T
}

namespace std {
    #define T Physica::Core::SIMD<Physica::Core::Real<Physica::Core::Float16>, 2>

    __host__ __device__ inline T max(T a, T b) noexcept {
        return T(__hmax2_nan(a.toMachine(), b.toMachine()));
    }

    __host__ __device__ inline T min(T a, T b) noexcept {
        return T(__hmin2_nan(a.toMachine(), b.toMachine()));
    }

    #undef T
}
