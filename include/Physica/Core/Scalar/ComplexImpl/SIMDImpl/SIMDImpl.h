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

#include "Physica/Core/Scalar/ComplexImpl/SIMD.h"

namespace Physica {
    template<Scalar T, int Size>
    SIMD<Complex<T>, Size>::SIMD(int x) : SIMD(ScalarType(x)) {}

    template<Scalar T, int Size>
    SIMD<Complex<T>, Size>::SIMD(double x) : SIMD(ScalarType(x)) {}

    template<Scalar T, int Size>
    SIMD<Complex<T>, Size>::SIMD(ScalarType x) {
        if constexpr (isSeparatable) {
            using Half = SIMD<Complex<T>, Size / 2>;
            const Half h(x);
            *this = asComplex({h.asReal(), h.asReal()});
        }
        else {
            if constexpr (T::Prec == Float32)
                *this = asComplex({x.real(), x.imag(), x.real(), x.imag()});
            else
                *this = asComplex({x.real(), x.imag()});
        }
    }

    template<Scalar T, int Size>
    SIMD<Complex<T>, Size>::SIMD(ScalarType x, int count) : SIMD(x) {
        cutoff(count);
    }

    template<Scalar T, int Size>
    SIMD<Complex<T>, Size>::ScalarType SIMD<Complex<T>, Size>::operator[](int index) const {
        return ScalarType(FullRealType::operator[](2 * index), FullRealType::operator[](2 * index + 1));
    }

    template<Scalar T, int Size>
    SIMD<Complex<T>, Size> SIMD<Complex<T>, Size>::operator+(const SIMD other) const {
        return asComplex(FullRealType::operator+(other));
    }

    template<Scalar T, int Size>
    SIMD<Complex<T>, Size> SIMD<Complex<T>, Size>::operator-(const SIMD other) const {
        return asComplex(FullRealType::operator-(other));
    }
    /**
     * Reference:
     * [1] vectorclass add-on; https://github.com/vectorclass/add-on
     */
    template<Scalar T, int Size>
    SIMD<Complex<T>, Size> SIMD<Complex<T>, Size>::operator*(const SIMD other) const {
        const auto pair = makeFullRealImag();
        return asComplex(mul_addsub(pair.first, other.asReal(), pair.second * other.swapRealImag()));
    }

    template<Scalar T, int Size>
    SIMD<Complex<T>, Size> SIMD<Complex<T>, Size>::operator*(const ScalarType& x) const {
        return operator*(SIMD(x));
    }

    template<Scalar T, int Size>
    SIMD<Complex<T>, Size> SIMD<Complex<T>, Size>::operator*(const T& x) const {
        return asComplex(FullRealType::operator*(x));
    }
    /**
     * Reference:
     * [1] vectorclass add-on; https://github.com/vectorclass/add-on
     */
    template<Scalar T, int Size>
    SIMD<Complex<T>, Size> SIMD<Complex<T>, Size>::operator/(const SIMD other) const {
        const auto pair = makeFullRealImag();
        const T factor = reciprocal(abs(other.asReal()).max()); // Avoid underflow
        const auto normed = other * factor;
        return asComplex(mul_addsub(pair.second, normed.swapRealImag(), -pair.first * normed.asReal()) / normed.squaredNorm() * factor);
    }

    template<Scalar T, int Size>
    SIMD<Complex<T>, Size> SIMD<Complex<T>, Size>::operator-() const {
        return asComplex(FullRealType::operator-());
    }

    template<Scalar T, int Size>
    void SIMD<Complex<T>, Size>::load(const ScalarType* p) {
        FullRealType::load(reinterpret_cast<const T*>(p));
    }

    template<Scalar T, int Size>
    void SIMD<Complex<T>, Size>::load_partial(const ScalarType* p, int n) {
        FullRealType::load_partial(reinterpret_cast<const T*>(p), 2 * n);
    }

    template<Scalar T, int Size>
    void SIMD<Complex<T>, Size>::store(ScalarType* p) const {
        FullRealType::store(reinterpret_cast<T*>(p));
    }

    template<Scalar T, int Size>
    void SIMD<Complex<T>, Size>::store_partial(ScalarType* p, int n) const {
        FullRealType::store_partial(reinterpret_cast<T*>(p), 2 * n);
    }

    template<Scalar T, int Size>
    auto SIMD<Complex<T>, Size>::cutoff(int count) -> SIMD& {
        assert(0 < count && count < int(Size) && "[Error]: Invalid count");
        RealBase::cutoff(2 * count);
        return *this;
    }

    template<Scalar T, int Size>
    auto SIMD<Complex<T>, Size>::conjugate() const noexcept -> This {
        if constexpr (Size == 1)
            return asComplex(asReal().template change_sign<0, 1>());
        else if constexpr (Size == 2)
            return asComplex(asReal().template change_sign<0, 1, 0, 1>());
        else if constexpr (Size == 4)
            return asComplex(asReal().template change_sign<0, 1, 0, 1, 0, 1, 0, 1>());
        else {
            static_assert(Size == 8, "[Error]: Unexpected type");
            return asComplex(asReal().template change_sign<0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1>());
        }
    }

    template<Scalar T, int Size>
    auto SIMD<Complex<T>, Size>::makeFullRealImag() const noexcept -> FullRealPair {
        FullRealType re;
        FullRealType im;
        if constexpr (T::Prec == Float32) {
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

    template<Scalar T, int Size>
    auto SIMD<Complex<T>, Size>::sum() const -> ScalarType {
        if constexpr (isSeparatable)
            return getHigh().sum() + getLow().sum();
        else {
            if constexpr (T::Prec == Float32)
                return operator[](0) + operator[](1);
            else
                return operator[](0);
        }
    }

    template<Scalar T, int Size>
    auto SIMD<Complex<T>, Size>::real() const noexcept -> RealType {
        if constexpr (isSeparatable)
            return permRealImag().getLow();
        else {
            const RealBase zero(0);
            if constexpr (T::Prec == Float32)
                return RealBase::template blend<0, 4, 2, 6>(asReal(), zero);
            else
                return RealBase::template blend<0, 2>(asReal(), zero);
        }
    }

    template<Scalar T, int Size>
    auto SIMD<Complex<T>, Size>::imag() const noexcept -> RealType {
        if constexpr (isSeparatable)
            return permRealImag().getHigh();
        else {
            const RealBase zero(0);
            if constexpr (T::Prec == Float32)
                return RealBase::template blend<1, 4, 3, 6>(asReal(), zero);
            else
                return RealBase::template blend<1, 3>(asReal(), zero);
        }
    }

    template<Scalar T, int Size>
    auto SIMD<Complex<T>, Size>::asComplex(FullRealType reals) -> SIMD {
        This result{};
        static_cast<FullRealType&>(result) = std::move(reals);
        return result;
    }

    template<Scalar T, int Size>
    SIMD<Complex<T>, Size> mul_add(const SIMD<Complex<T>, Size> a, const SIMD<Complex<T>, Size> b, const SIMD<Complex<T>, Size> c) noexcept {
        return a * b + c;
    }
}
