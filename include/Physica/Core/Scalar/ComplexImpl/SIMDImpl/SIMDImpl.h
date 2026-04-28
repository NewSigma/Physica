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
        return ScalarType(storage[2 * index], storage[2 * index + 1]);
    }

    template<Scalar T, int Size>
    auto SIMD<Complex<T>, Size>::operator+(const Packet auto x) const {
        static_assert(x.size() == size(), "[Error]: Size mismatch");
        if constexpr (x.isComplex())
            return asComplex(asReal() + x.asReal());
        else
            return asComplex(asReal() + FullRealType(x, SIMD<T, Size>(0)).gatherRealImag());
    }

    template<Scalar T, int Size>
    auto SIMD<Complex<T>, Size>::operator-(const Packet auto x) const {
        static_assert(x.size() == size(), "[Error]: Size mismatch");
        if constexpr (x.isComplex())
            return asComplex(asReal() - x.asReal());
        else
            return asComplex(asReal() - FullRealType(x, SIMD<T, Size>(0)).gatherRealImag());
    }

    template<Scalar T, int Size>
    auto SIMD<Complex<T>, Size>::operator*(Packet auto x) const {
        static_assert(x.size() == size(), "[Error]: Size mismatch");
        if constexpr (x.isComplex()) {
            /**
             * Reference:
             * [1] vectorclass add-on; https://github.com/vectorclass/add-on
             */
            const auto pair = makeFullRealImag();
            return asComplex(mul_addsub(pair.first, x.asReal(), pair.second * x.swapRealImag()));
        }
        else
            return asComplex(asReal() * FullRealType(x, x).gatherRealImag());
    }

    template<Scalar T, int Size>
    auto SIMD<Complex<T>, Size>::operator*(const Scalar auto& x) const {
        if constexpr (x.isComplex())
            return operator*(SIMD(x));
        else
            return asComplex(storage * x);
    }
    /**
     * Reference:
     * [1] vectorclass add-on; https://github.com/vectorclass/add-on
     */
    template<Scalar T, int Size>
    auto SIMD<Complex<T>, Size>::operator/(const Packet auto x) const {
        static_assert(x.size() == size(), "[Error]: Size mismatch");
        if constexpr (x.isComplex()) {
            const auto pair = makeFullRealImag();
            const T factor = reciprocal(abs(x.asReal()).max()); // Avoid underflow
            const auto normed = x * factor;
            return asComplex(mul_addsub(pair.second, normed.swapRealImag(), -pair.first * normed.asReal()) / normed.squaredNorm() * factor);
        }
        else
            return asComplex(asReal() / FullRealType(x, x).gatherRealImag());
    }

    template<Scalar T, int Size>
    SIMD<Complex<T>, Size> SIMD<Complex<T>, Size>::operator-() const {
        return asComplex(-storage);
    }

    template<Scalar T, int Size>
    void SIMD<Complex<T>, Size>::load(const ScalarType* p) noexcept {
        storage.load(reinterpret_cast<const T*>(p));
    }

    template<Scalar T, int Size>
    void SIMD<Complex<T>, Size>::load(const ScalarType* p, int n) noexcept {
        assert(0 < n && n < Size && "[Error]: Invalid size for partial operation");
        storage.load(reinterpret_cast<const T*>(p), 2 * n);
    }

    template<Scalar T, int Size>
    void SIMD<Complex<T>, Size>::store(ScalarType* p) const noexcept {
        storage.store(reinterpret_cast<T*>(p));
    }

    template<Scalar T, int Size>
    void SIMD<Complex<T>, Size>::store(ScalarType* p, int n) const noexcept {
        assert(0 < n && n < Size && "[Error]: Invalid size for partial operation");
        storage.store(reinterpret_cast<T*>(p), 2 * n);
    }

    template<Scalar T, int Size>
    auto SIMD<Complex<T>, Size>::cutoff(int count) -> SIMD& {
        assert(0 < count && count < int(Size) && "[Error]: Invalid count");
        storage.cutoff(2 * count);
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
            re = storage.template shuffle<0, 0, 2, 2>();
            im = storage.template shuffle<1, 1, 3, 3>();
        }
        else {
            if constexpr (Size == 1) {
                re = storage.template shuffle<0, 0>();
                im = storage.template shuffle<1, 1>();
            }
            else if constexpr (Size == 2) {
                re = storage.template shuffle<0, 0, 0, 0>();
                im = storage.template shuffle<1, 1, 1, 1>();
            }
            else {
                static_assert(Size == 4, "[Error]: Unexpected size");
                re = storage.template shuffle<0, 0, 0, 0, 0, 0, 0, 0>();
                im = storage.template shuffle<1, 1, 1, 1, 1, 1, 1, 1>();
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
    auto SIMD<Complex<T>, Size>::isZero() const noexcept {
        return storage.isZero();
    }

    template<Scalar T, int Size>
    auto SIMD<Complex<T>, Size>::isFinite() const noexcept {
        return storage.isFinite();
    }

    template<Scalar T, int Size>
    auto SIMD<Complex<T>, Size>::isSubNormal() const noexcept {
        auto bools = Base::scatterRealImag().isSubNormal();
        return bools.getLow() && bools.getHigh();
    }

    template<Scalar T, int Size>
    auto SIMD<Complex<T>, Size>::real() const noexcept -> RealType {
        if constexpr (isSeparatable)
            return gatherRealImag().getLow();
        else if constexpr (T::Prec == Float32)
            return storage.template blend<0, 4, 2, 6>(asReal(), FullRealType::zeros());
        else
            return storage.template blend<0, 2>(asReal(), FullRealType::zeros());
    }

    template<Scalar T, int Size>
    auto SIMD<Complex<T>, Size>::imag() const noexcept -> RealType {
        if constexpr (isSeparatable)
            return gatherRealImag().getHigh();
        else if constexpr (T::Prec == Float32)
            return storage.template blend<1, 4, 3, 6>(asReal(), FullRealType::zeros());
        else
            return storage.template blend<1, 3>(asReal(), FullRealType::zeros());
    }

    template<Scalar T, int Size>
    auto SIMD<Complex<T>, Size>::zeros() noexcept -> SIMD {
        return asComplex(FullRealType::zeros());
    }

    template<Scalar T, int Size>
    auto SIMD<Complex<T>, Size>::asComplex(FullRealType storage) noexcept -> SIMD {
        SIMD result{};
        result.storage = storage;
        return result;
    }

    template<Scalar T, Scalar U, int Size>
    __host__ __device__ auto operator/(const SIMD<T, Size> x, const SIMD<Complex<U>, Size> y) requires(!T::isComplex() && !T::isDiffable()) {
        return y.conjugate() * x / SIMD<Complex<U>, Size>::asComplex(y.squaredNorm()).real();
    }
}
