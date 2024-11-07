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
    template<class T, size_t Size>
    SIMD<Complex<T>, Size>::SIMD(int x) : SIMD(ScalarType(x)) {}

    template<class T, size_t Size>
    SIMD<Complex<T>, Size>::SIMD(ScalarType x) {
        if constexpr (isSeparatable) {
            using Half = SIMD<Complex<T>, Size / 2>;
            const Half h(x);
            *this = This(Base(h.getImpl(), h.getImpl()));
        }
        else {
            if constexpr (T::Option == Float32)
                *this = This(Base(x.real().toMachine(), x.imag().toMachine(), x.real().toMachine(), x.imag().toMachine()));
            else
                *this = This(Base(x.real().toMachine(), x.imag().toMachine()));
        }
    }

    template<class T, size_t Size>
    SIMD<Complex<T>, Size>::SIMD(Base base) : Base(std::move(base)) {}

    template<class T, size_t Size>
    inline typename SIMD<Complex<T>, Size>::ScalarType SIMD<Complex<T>, Size>::operator[](int index) const {
        return ScalarType(Base::operator[](2 * index), Base::operator[](2 * index + 1));
    }

    template<class T, size_t Size>
    inline SIMD<Complex<T>, Size> SIMD<Complex<T>, Size>::operator+(const SIMD& other) const {
        return This(Base::operator+(other));
    }

    template<class T, size_t Size>
    inline SIMD<Complex<T>, Size> SIMD<Complex<T>, Size>::operator-(const SIMD& other) const {
        return This(Base::operator-(other));
    }
    /**
     * Reference:
     * [1] vectorclass add-on; https://github.com/vectorclass/add-on
     */
    template<class T, size_t Size>
    inline SIMD<Complex<T>, Size> SIMD<Complex<T>, Size>::operator*(const SIMD& other) const {
        Base a_re, a_im;
        if constexpr (T::Option == Float32) {
            a_re = Base::template shuffle<0, 0, 2, 2>();
            a_im = Base::template shuffle<1, 1, 3, 3>();
        }
        else {
            if constexpr (Size == 1) {
                a_re = Base::template shuffle<0, 0>();
                a_im = Base::template shuffle<1, 1>();
            }
            else if constexpr (Size == 2) {
                a_re = Base::template shuffle<0, 0, 0, 0>();
                a_im = Base::template shuffle<1, 1, 1, 1>();
            }
            else {
                static_assert(Size == 4, "[Error]: Unexpected size");
                a_re = Base::template shuffle<0, 0, 0, 0, 0, 0, 0, 0>();
                a_im = Base::template shuffle<1, 1, 1, 1, 1, 1, 1, 1>();
            }
        }
        return This(mul_addsub(a_re, other.getImpl(), a_im * other.swapRealImag()));
    }

    template<class T, size_t Size>
    inline SIMD<Complex<T>, Size> SIMD<Complex<T>, Size>::operator*(const ScalarType& x) const {
        return operator*(SIMD(x));
    }

    template<class T, size_t Size>
    inline SIMD<Complex<T>, Size> SIMD<Complex<T>, Size>::operator*(const T& x) const {
        return This(Base::operator*(x));
    }

    template<class T, size_t Size>
    inline SIMD<Complex<T>, Size> SIMD<Complex<T>, Size>::operator-() const {
        return This(Base::operator-());
    }

    template<class T, size_t Size>
    inline void SIMD<Complex<T>, Size>::load(const ScalarType* p) {
        Base::load(reinterpret_cast<const T*>(p));
    }

    template<class T, size_t Size>
    inline void SIMD<Complex<T>, Size>::load_partial(const ScalarType* p, int n) {
        Base::load_partial(reinterpret_cast<const T*>(p), 2 * n);
    }

    template<class T, size_t Size>
    inline void SIMD<Complex<T>, Size>::store(ScalarType* p) const {
        Base::store(reinterpret_cast<T*>(p));
    }

    template<class T, size_t Size>
    inline void SIMD<Complex<T>, Size>::store_partial(ScalarType* p, int n) const {
        Base::store_partial(reinterpret_cast<T*>(p), 2 * n);
    }

    template<class T, size_t Size>
    typename SIMD<Complex<T>, Size>::Base SIMD<Complex<T>, Size>::swapRealImag() const noexcept {
        Base result;
        if constexpr (T::Option == Float32)
            result = Base::template shuffle<1, 0, 3, 2>();
        else {
            if constexpr (Size == 1)
                result = Base::template shuffle<1, 0>();
            else if constexpr (Size == 2)
                result = Base::template shuffle<1, 0, 1, 0>();
            else {
                static_assert(Size == 4, "[Error]: Unexpected size");
                result = Base::template shuffle<1, 0, 1, 0, 1, 0, 1, 0>();
            }
        }
        return result;
    }

    template<class T, size_t Size>
    inline typename SIMD<Complex<T>, Size>::Base SIMD<Complex<T>, Size>::permRealImag() const noexcept {
        auto result = getImpl();
        if constexpr (Size == 2)
            result = result.template permute<0, 2, 1, 3>();
        else if constexpr (Size == 4)
            result = result.template permute<0, 2, 4, 6, 1, 3, 5, 7>();
        else if constexpr (Size == 8)
            result = result.template permute<0, 2, 4, 6, 8, 10, 12, 14, 1, 3, 5, 7, 9, 11, 13, 15>();
        else
            static_assert(Size == 1, "[Error]: Unexpected size");
        return result;
    }

    template<class T, size_t Size>
    inline typename SIMD<Complex<T>, Size>::ScalarType SIMD<Complex<T>, Size>::sum() const {
        if constexpr (isSeparatable)
            return getHigh().sum() + getLow().sum();
        else {
            if constexpr (T::Option == Float32)
                return operator[](0) + operator[](1);
            else
                return operator[](0);
        }
    }

    template<class T, size_t Size>
    inline typename SIMD<Complex<T>, Size>::RealType SIMD<Complex<T>, Size>::real() const noexcept {
        return permRealImag().getLow();
    }

    template<class T, size_t Size>
    inline typename SIMD<Complex<T>, Size>::RealType SIMD<Complex<T>, Size>::imag() const noexcept {
        return permRealImag().getHigh();
    }
}
