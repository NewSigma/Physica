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
    template<Scalar T, size_t Size>
    SIMD<Complex<T>, Size>::SIMD(int x) : SIMD(ScalarType(x)) {}

    template<Scalar T, size_t Size>
    SIMD<Complex<T>, Size>::SIMD(ScalarType x) {
        if constexpr (isSeparatable) {
            using Half = SIMD<Complex<T>, Size / 2>;
            const Half h(x);
            *this = asComplex({h.asReal(), h.asReal()});
        }
        else {
            if constexpr (T::Option == Float32)
                *this = asComplex({x.real().toMachine(), x.imag().toMachine(), x.real().toMachine(), x.imag().toMachine()});
            else
                *this = asComplex({x.real().toMachine(), x.imag().toMachine()});
        }
    }

    template<Scalar T, size_t Size>
    inline SIMD<Complex<T>, Size>::ScalarType SIMD<Complex<T>, Size>::operator[](int index) const {
        return ScalarType(FullRealType::operator[](2 * index), FullRealType::operator[](2 * index + 1));
    }

    template<Scalar T, size_t Size>
    inline SIMD<Complex<T>, Size> SIMD<Complex<T>, Size>::operator+(const SIMD& other) const {
        return asComplex(FullRealType::operator+(other));
    }

    template<Scalar T, size_t Size>
    inline SIMD<Complex<T>, Size> SIMD<Complex<T>, Size>::operator-(const SIMD& other) const {
        return asComplex(FullRealType::operator-(other));
    }
    /**
     * Reference:
     * [1] vectorclass add-on; https://github.com/vectorclass/add-on
     */
    template<Scalar T, size_t Size>
    inline SIMD<Complex<T>, Size> SIMD<Complex<T>, Size>::operator*(const SIMD& other) const {
        const auto pair = makeFullReIm();
        return asComplex(mul_addsub(pair.first, other.asReal(), pair.second * other.swapRealImag()));
    }

    template<Scalar T, size_t Size>
    inline SIMD<Complex<T>, Size> SIMD<Complex<T>, Size>::operator*(const ScalarType& x) const {
        return operator*(SIMD(x));
    }

    template<Scalar T, size_t Size>
    inline SIMD<Complex<T>, Size> SIMD<Complex<T>, Size>::operator*(const T& x) const {
        return asComplex(FullRealType::operator*(x));
    }
    /**
     * Reference:
     * [1] vectorclass add-on; https://github.com/vectorclass/add-on
     */
    template<Scalar T, size_t Size>
    inline SIMD<Complex<T>, Size> SIMD<Complex<T>, Size>::operator/(const SIMD& other) const {
        const auto pair = makeFullReIm();
        return asComplex(mul_addsub(pair.second, other.swapRealImag(), -pair.first * other.asReal()) / other.squaredNorm());
    }

    template<Scalar T, size_t Size>
    inline SIMD<Complex<T>, Size> SIMD<Complex<T>, Size>::operator-() const {
        return asComplex(FullRealType::operator-());
    }

    template<Scalar T, size_t Size>
    inline void SIMD<Complex<T>, Size>::load(const ScalarType* p) {
        FullRealType::load(reinterpret_cast<const T*>(p));
    }

    template<Scalar T, size_t Size>
    inline void SIMD<Complex<T>, Size>::load_partial(const ScalarType* p, int n) {
        FullRealType::load_partial(reinterpret_cast<const T*>(p), 2 * n);
    }

    template<Scalar T, size_t Size>
    inline void SIMD<Complex<T>, Size>::store(ScalarType* p) const {
        FullRealType::store(reinterpret_cast<T*>(p));
    }

    template<Scalar T, size_t Size>
    inline void SIMD<Complex<T>, Size>::store_partial(ScalarType* p, int n) const {
        FullRealType::store_partial(reinterpret_cast<T*>(p), 2 * n);
    }

    template<Scalar T, size_t Size>
    inline SIMD<Complex<T>, Size>::FullRealType SIMD<Complex<T>, Size>::permRealImag() const noexcept {
        auto result = asReal();
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

    template<Scalar T, size_t Size>
    inline SIMD<Complex<T>, Size>::ScalarType SIMD<Complex<T>, Size>::sum() const {
        if constexpr (isSeparatable)
            return getHigh().sum() + getLow().sum();
        else {
            if constexpr (T::Option == Float32)
                return operator[](0) + operator[](1);
            else
                return operator[](0);
        }
    }

    template<Scalar T, size_t Size>
    inline SIMD<Complex<T>, Size>::RealType SIMD<Complex<T>, Size>::real() const noexcept {
        return permRealImag().getLow();
    }

    template<Scalar T, size_t Size>
    inline SIMD<Complex<T>, Size>::RealType SIMD<Complex<T>, Size>::imag() const noexcept {
        return permRealImag().getHigh();
    }

    template<Scalar T, size_t Size>
    SIMD<Complex<T>, Size> SIMD<Complex<T>, Size>::asComplex(FullRealType reals) {
        This result{};
        static_cast<FullRealType&>(result) = std::move(reals);
        return result;
    }

    template<Scalar T, size_t Size>
    inline SIMD<Complex<T>, Size>::FullRealPair SIMD<Complex<T>, Size>::makeFullReIm() const noexcept {
        FullRealType re, im;
        if constexpr (T::Option == Float32) {
            re = FullRealType::template shuffle<0, 0, 2, 2>();
            im = FullRealType::template shuffle<1, 1, 3, 3>();
        }
        else {
            if constexpr (Size == 1) {
                re = FullRealType::template shuffle<0, 0>();
                im = FullRealType::template shuffle<1, 1>();
            }
            else if constexpr (Size == 2) {
                re = FullRealType::template shuffle<0, 0, 0, 0>();
                im = FullRealType::template shuffle<1, 1, 1, 1>();
            }
            else {
                static_assert(Size == 4, "[Error]: Unexpected size");
                re = FullRealType::template shuffle<0, 0, 0, 0, 0, 0, 0, 0>();
                im = FullRealType::template shuffle<1, 1, 1, 1, 1, 1, 1, 1>();
            }
        }
        return std::make_pair(std::move(re), std::move(im));
    }
}
