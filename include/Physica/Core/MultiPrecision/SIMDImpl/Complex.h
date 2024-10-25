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
        using RealPacket = device_obj<BestPacket<ScalarType, Length * 2>>;
    public:
        constexpr static size_t Size = 1;
        using Type = typename RealPacket::Type;
    };

    template<class T, size_t Size>
    class SIMD<Complex<T>, Size> : private SIMD<T, Size * 2> {
        static_assert(T::Option == Float32 || T::Option == Float64, "[Error]: Not implemented");
        using ScalarType = Complex<T>;
        using This = SIMD<ScalarType, Size>;
        using Base = SIMD<T, Size * 2>;
        using BoolSIMDType = typename Base::BoolSIMDType;
        using TrivialType = typename ScalarType::TrivialType;
        using HalfType = typename std::conditional<sizeof(Base) * CHAR_BIT != 128, SIMD<Complex<T>, Size / 2>, PlainStruct<void>>::type;
        using Base::isSeparatable;
    public:
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
        [[nodiscard]] inline SIMD operator+(const SIMD& other) const { return Base::operator+(other); }
        [[nodiscard]] inline SIMD operator-(const SIMD& other) const { return Base::operator-(other); }
        [[nodiscard]] inline SIMD operator*(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator*(const ScalarType& x) const { return operator*(SIMD(x)); }
        [[nodiscard]] inline SIMD operator*(const T& x) const { return Base::operator*(x); }
        //[[nodiscard]] inline SIMD operator/(const SIMD& other) const;
        [[nodiscard]] inline SIMD operator-() const { return Base::operator-(); }
        void operator+=(const SIMD& other) { *this = *this + other; }
        void operator-=(const SIMD& other) { *this = *this - other; }
        void operator*=(const SIMD& other) { *this = *this * other; }
        //void operator/=(const SIMD& other) { *this = *this / other; }
        /* Operations */
        inline void load(const ScalarType* p);
        inline void load_partial(int n, const ScalarType* p);
        inline void store(ScalarType* p) const;
        inline void store_partial(int n, ScalarType* p) const;
        //inline void insert(int index, const ScalarType& value);
        [[nodiscard]] inline ScalarType horizontal_add() const;
        void swap(SIMD& __restrict other) noexcept { std::swap(*this, other); }
        /* Getters */
        [[nodiscard]] constexpr static size_t size() { return Size; }
        [[nodiscard]] Base& getImpl() noexcept { return *this; }
        [[nodiscard]] const Base& getImpl() const noexcept { return *this; }
        [[nodiscard]] HalfType getLow() const noexcept { return HalfType::asComplex(Base::getLow()); }
        [[nodiscard]] HalfType getHigh() const noexcept { return HalfType::asComplex(Base::getHigh()); }
        /* Static members */
        [[nodiscard]] static SIMD asComplex(Base base) noexcept { return base; }
        template<class RandomType>
        [[nodiscard]] static SIMD random_uniform(RandomType& gen) { return asComplex(Base::random_uniform(gen)); }
    private:
        SIMD(Base base) : Base(std::move(base)) {}
    };

    template<class T, size_t Size>
    SIMD<Complex<T>, Size>::SIMD(int x) : SIMD(ScalarType(x)) {}

    template<class T, size_t Size>
    SIMD<Complex<T>, Size>::SIMD(ScalarType x) {
        if constexpr (isSeparatable) {
            using Half = SIMD<Complex<T>, Size / 2>;
            const Half h(x);
            *this = Base(h.getImpl(), h.getImpl());
        }
        else {
            if constexpr (T::Option == Float32)
                *this = Base(x.real().getTrivial(), x.imag().getTrivial(), x.real().getTrivial(), x.imag().getTrivial());
            else
                *this = Base(x.real().getTrivial(), x.imag().getTrivial());
        }
    }

    template<class T, size_t Size>
    inline typename SIMD<Complex<T>, Size>::ScalarType SIMD<Complex<T>, Size>::operator[](int index) const {
        return ScalarType(Base::operator[](2 * index), Base::operator[](2 * index + 1));
    }
    /**
     * Reference:
     * [1] vectorclass add-on; https://github.com/vectorclass/add-on
     */
    template<class T, size_t Size>
    inline SIMD<Complex<T>, Size> SIMD<Complex<T>, Size>::operator*(const SIMD& other) const {
        Base a_re, a_im, b_flip;
        if constexpr (T::Option == Float32) {
            a_re = Base::template shuffle<0, 0, 2, 2>();
            a_im = Base::template shuffle<1, 1, 3, 3>();
            b_flip = other.template shuffle<1, 0, 3, 2>();
        }
        else {
            if constexpr (Size == 1) {
                a_re = Base::template shuffle<0, 0>();
                a_im = Base::template shuffle<1, 1>();
                b_flip = other.template shuffle<1, 0>();
            }
            else if constexpr (Size == 2) {
                a_re = Base::template shuffle<0, 0, 0, 0>();
                a_im = Base::template shuffle<1, 1, 1, 1>();
                b_flip = other.template shuffle<1, 0, 1, 0>();
            }
            else {
                static_assert(Size == 4, "[Error]: Unexpected size");
                a_re = Base::template shuffle<0, 0, 0, 0, 0, 0, 0, 0>();
                a_im = Base::template shuffle<1, 1, 1, 1, 1, 1, 1, 1>();
                b_flip = other.template shuffle<1, 0, 1, 0, 1, 0, 1, 0>();
            }
        }
        return mul_addsub(a_re, other.getImpl(), a_im * b_flip);
    }

    template<class T, size_t Size>
    inline void SIMD<Complex<T>, Size>::load(const ScalarType* p) {
        Base::load(reinterpret_cast<const T*>(p));
    }

    template<class T, size_t Size>
    inline void SIMD<Complex<T>, Size>::load_partial(int n, const ScalarType* p) {
        Base::load_partial(2 * n, reinterpret_cast<const T*>(p));
    }

    template<class T, size_t Size>
    inline void SIMD<Complex<T>, Size>::store(ScalarType* p) const {
        Base::store(reinterpret_cast<T*>(p));
    }

    template<class T, size_t Size>
    inline void SIMD<Complex<T>, Size>::store_partial(int n, ScalarType* p) const {
        Base::store_partial(2 * n, reinterpret_cast<T*>(p));
    }

    template<class T, size_t Size>
    inline typename SIMD<Complex<T>, Size>::ScalarType SIMD<Complex<T>, Size>::horizontal_add() const {
        if constexpr (isSeparatable)
            return getHigh().horizontal_add() + getLow().horizontal_add();
        else {
            if constexpr (T::Option == Float32)
                return operator[](0) + operator[](1);
            else
                return operator[](0);
        }
    }
}

namespace Physica {
    template<class T, size_t Size>
    class Traits<Core::SIMD<Core::Complex<T>, Size>> {
    public:
        using ScalarType = Core::Complex<T>;
    };
}

#include "MathComplex.h"
