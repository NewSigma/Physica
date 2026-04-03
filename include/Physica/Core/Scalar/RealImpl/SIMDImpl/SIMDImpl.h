/*
 * Copyright 2023-2026 Weibo He.
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

#include "Instruset.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Scalar/RealImpl/SIMD.h"
#include "Physica/Core/Utils/Builtin.h"

namespace Physica {
    template<Scalar T, int Size>
    constexpr SIMD<T, Size>::SIMD(T::MachineType x) noexcept : pack(x) {}

    template<Scalar T, int Size>
    constexpr SIMD<T, Size>::SIMD(T x) noexcept : pack(x.toMachine()) {}

    template<Scalar T, int Size>
    SIMD<T, Size>::SIMD(T x, int count) noexcept {
        assert(0 < count && count <= int(Size) && "[Error]: Invalid count");
        if (count == Size) {
            *this = SIMD(x);
            return;
        }

        constexpr int HalfSize = Size / 2;
        if constexpr (isSeparatable) {
            if (count > HalfSize)
                *this = SIMD(HalfType(x), HalfType(x, count - HalfSize));
            else
                *this = SIMD(HalfType(x, count), HalfType(0));
        }
        else {
            if constexpr (T::Prec == Float32) {
                const float f = float(x);
                switch (count) {
                case 1:
                    *this = SIMD(Pack(f, 0.0F, 0.0F, 0.0F));
                    break;
                case 2:
                    *this = SIMD(Pack(f, f, 0.0F, 0.0F));
                    break;
                case 3:
                    *this = SIMD(Pack(f, f, f, 0.0F));
                    break;
                default:
                    unreachable();
                }
            }
            else {
                static_assert(T::Prec == Float64, "[Error]: Unsupported float type");
                *this = SIMD(Pack(double(x), 0.0));
            }
        }
    }

    template<Scalar T, int Size>
    SIMD<T, Size>::SIMD(Scalar auto... args) noexcept : pack(args.toMachine()...) {
        static_assert(sizeof...(args) == Size, "[Error]: Number of elements does not match");
    }

    template<Scalar T, int Size>
    SIMD<T, Size>::SIMD(HalfType a, HalfType b) noexcept : pack(a.toMachine(), b.toMachine()) {}

    template<Scalar T, int Size>
    SIMD<T, Size>::SIMD(BoolSIMDType x) noexcept : pack(x.toMachine()) {}

    template<Scalar T, int Size>
    T SIMD<T, Size>::operator[](int index) const {
        if constexpr (isForward)
            return T(pack.operator[](index * 2), pack.operator[](index * 2 + 1));
        else
            return T(pack.operator[](index));
    }

    template<Scalar T, int Size>
    SIMD<T, Size> SIMD<T, Size>::operator+(const SIMD other) const {
        return SIMD(toMachine() + other.toMachine());
    }

    template<Scalar T, int Size>
    SIMD<T, Size> SIMD<T, Size>::operator-(const SIMD other) const {
        return SIMD(toMachine() - other.toMachine());
    }

    template<Scalar T, int Size>
    SIMD<T, Size> SIMD<T, Size>::operator*(const SIMD other) const {
        return SIMD(toMachine() * other.toMachine());
    }

    template<Scalar T, int Size>
    SIMD<T, Size> SIMD<T, Size>::operator*(const T& x) const {
        return operator*(SIMD(x));
    }

    template<Scalar T, int Size>
    SIMD<T, Size> SIMD<T, Size>::operator/(const SIMD other) const {
        // Vector partial reciprocal might divide by zero, do not assert it.
        return SIMD(toMachine() / other.toMachine());
    }

    template<Scalar T, int Size>
    SIMD<T, Size> SIMD<T, Size>::operator-() const {
        return SIMD(-toMachine());
    }

    template<Scalar T, int Size>
    SIMD<T, Size> SIMD<T, Size>::operator&(const SIMD other) const {
        return SIMD(toMachine() & other.toMachine());
    }

    template<Scalar T, int Size>
    SIMD<T, Size> SIMD<T, Size>::operator|(const SIMD other) const {
        return SIMD(toMachine() | other.toMachine());
    }

    template<Scalar T, int Size>
    SIMD<T, Size> SIMD<T, Size>::operator^(const SIMD other) const {
        return SIMD(toMachine() ^ other.toMachine());
    }

    template<Scalar T, int Size>
    SIMD<T, Size> SIMD<T, Size>::operator!() const {
        return SIMD(!toMachine());
    }

    template<Scalar T, int Size>
    auto SIMD<T, Size>::operator==(const SIMD other) const {
        return BoolSIMDType(toMachine() == other.toMachine());
    }

    template<Scalar T, int Size>
    auto SIMD<T, Size>::operator>(const SIMD other) const {
        return BoolSIMDType(toMachine() > other.toMachine());
    }

    template<Scalar T, int Size>
    auto SIMD<T, Size>::operator<(const SIMD other) const {
        return BoolSIMDType(toMachine() < other.toMachine());
    }

    template<Scalar T, int Size>
    void SIMD<T, Size>::load(const T* p) & noexcept {
        pack.load(reinterpret_cast<const typename T::MachineType*>(p));
    }
    /**
     * Adapted from [1], we require a more strict range of \param n.
     * SIMD<T, Size>::store() is similar
     *
     * Reference:
     * [1] vectorclass2; https://github.com/vectorclass/version2
     */
    template<Scalar T, int Size>
    void SIMD<T, Size>::load(const T* p, int n) & noexcept {
        assert(0 < n && n < Size && "[Error]: Invalid size for partial operation");
        if constexpr (T::Prec == Float32) {
            if constexpr (Size == 4) {
                if constexpr (Instruset::hasAVX512VL())
                    *this = _mm_maskz_loadu_ps(__mmask8((1U << n) - 1), (float*)p);
                else {
                    switch (n) {
                    case 1:
                        *this = _mm_load_ss((float*)p);
                        break;
                    case 2:
                        *this = _mm_castpd_ps(_mm_load_sd((double*)p));
                        break;
                    case 3: {
                        auto t1 = _mm_castpd_ps(_mm_load_sd((double*)p));
                        auto t2 = _mm_load_ss((float*)p + 2);
                        *this = _mm_movelh_ps(t1, t2);
                        break;
                    }
                    default:
                        unreachable();
                    }
                }
            }
            else if constexpr (Size == 8) {
                if constexpr (Instruset::hasAVX512VL())
                    *this = _mm256_maskz_loadu_ps(__mmask8((1U << n) - 1), (float*)p);
                else {
                    HalfType low{}, high{};
                    if (n < 4) {
                        low.load(p, n);
                        high = HalfType::zeros();
                    }
                    else if (n == 4) {
                        low.load(p);
                        high = HalfType::zeros();
                    }
                    else {
                        low.load(p);
                        high.load(p + 4, n - 4);
                    }
                    *this = This(low, high);
                }
            }
            else {
                static_assert(Size == 16, "[Error]: Unexpected size");
                *this = _mm512_maskz_loadu_ps(__mmask16((1 << n) - 1), (float*)p);
            }
        }
        else {
            static_assert(T::Prec == Float64);
            if constexpr (Size == 2)
                *this = _mm_load_sd((double*)p);
            else if constexpr (Size == 4) {
                if constexpr (Instruset::hasAVX512VL())
                    *this = _mm256_maskz_loadu_pd(__mmask8((1U << n) - 1), (double*)p);
                else {
                    HalfType low{}, high{};
                    if (n < 2) {
                        low.load(p, 1);
                        high = HalfType::zeros();
                    }
                    else if (n == 2) {
                        low.load(p);
                        high = HalfType::zeros();
                    }
                    else {
                        low.load(p);
                        high.load(p + 2, 1);
                    }
                    *this = This(low, high);
                }
            }
            else {
                static_assert(Size == 8, "[Error]: Unexpected size");
                *this = _mm512_maskz_loadu_pd(__mmask16((1 << n) - 1), (double*)p);
            }
        }
    }

    template<Scalar T, int Size>
    void SIMD<T, Size>::store(T* p) const noexcept {
        pack.store(reinterpret_cast<typename T::MachineType*>(p));
    }

    template<Scalar T, int Size>
    void SIMD<T, Size>::store(T* p, int n) const noexcept {
        assert(0 < n && n < Size && "[Error]: Invalid size for partial operation");
        if constexpr (T::Prec == Float32) {
            if constexpr (Size == 4) {
                if constexpr (Instruset::hasAVX512VL())
                    _mm_mask_storeu_ps((float*)p, __mmask8((1U << n) - 1), *this);
                else {
                    switch (n) {
                    case 1:
                        _mm_store_ss((float*)p, *this);
                        break;
                    case 2:
                        _mm_store_sd((double*)p, _mm_castps_pd(*this));
                        break;
                    case 3: {
                        _mm_store_sd((double*)p, _mm_castps_pd(*this));
                        _mm_store_ss((float*)p + 2, _mm_movehl_ps(*this, *this));
                        break;
                    }
                    default:
                        unreachable();
                    }
                }
            }
            else if constexpr (Size == 8) {
                if constexpr (Instruset::hasAVX512VL())
                    _mm256_mask_storeu_ps((float*)p, __mmask8((1U << n) - 1), *this);
                else {
                    if (n < 4)
                        getLow().store(p, n);
                    else if (n == 4)
                        getLow().store(p);
                    else {
                        getLow().store(p);
                        getHigh().store(p + 4, n - 4);
                    }
                }
            }
            else {
                static_assert(Size == 16, "[Error]: Unexpected size");
                _mm512_mask_storeu_ps((float*)p, __mmask16((1 << n) - 1), *this);
            }
        }
        else {
            static_assert(T::Prec == Float64);
            if constexpr (Size == 2)
                _mm_store_sd((double*)p, *this);
            else if constexpr (Size == 4) {
                if constexpr (Instruset::hasAVX512VL())
                    _mm256_mask_storeu_pd((double*)p, __mmask8((1U << n) - 1), *this);
                else {
                    if (n < 2)
                        getLow().store(p, 1);
                    else if (n == 2)
                        getLow().store(p);
                    else {
                        getLow().store(p);
                        getHigh().store(p + 2, 1);
                    }
                }
            }
            else {
                static_assert(Size == 8, "[Error]: Unexpected size");
                _mm512_mask_storeu_pd((double*)p, __mmask16((1 << n) - 1), *this);
            }
        }
    }

    template<Scalar T, int Size>
    void SIMD<T, Size>::insert(int index, const T& value) {
        if constexpr (isForward) {
            pack.insert(index * 2, value.value().toMachine());
            pack.insert(index * 2 + 1, value.grad().toMachine());
        }
        else
            pack.insert(index, value.toMachine());
    }

    template<Scalar T, int Size>
    template<int... Order>
    SIMD<T, Size> SIMD<T, Size>::shuffle() const {
        constexpr unsigned int mask = makeShuffleMask<Order...>(0);
        if constexpr (T::Prec == Float32) {
            static_assert(sizeof...(Order) == 4, "[Error]: Invalid number of orders");
            if constexpr (Size == 4)
                return _mm_shuffle_ps(*this, *this, mask);
            else if constexpr (Size == 8)
                return _mm256_shuffle_ps(*this, *this, mask);
            else {
                static_assert(Size == 16, "[Error]: Unexpected size");
                return _mm512_shuffle_ps(*this, *this, mask);
            }
        }
        else {
            static_assert(sizeof...(Order) == Size, "[Error]: Invalid number of orders");
            static_assert(T::Prec == Float64);
            if constexpr (Size == 2)
                return _mm_shuffle_pd(*this, *this, mask);
            else if constexpr (Size == 4)
                return _mm256_shuffle_pd(*this, *this, mask);
            else {
                static_assert(Size == 8, "[Error]: Unexpected size");
                return _mm512_shuffle_pd(*this, *this, mask);
            }
        }
    }

    template<Scalar T, int Size>
    template<int... Order>
    SIMD<T, Size> SIMD<T, Size>::permute() const {
        static_assert(sizeof...(Order) == Size, "[Error]: Size of Order do not match the packet");
        if constexpr (Size == 2)
            return Physica::permute2<Order...>(toMachine());
        else if constexpr (Size == 4)
            return Physica::permute4<Order...>(toMachine());
        else if constexpr (Size == 8)
            return Physica::permute8<Order...>(toMachine());
        else {
            static_assert(Size == 16, "[Error]: Unexpected size");
            return Physica::permute16<Order...>(toMachine());
        }
    }

    template<Scalar T, int Size>
    template<int... Flags>
    SIMD<T, Size> SIMD<T, Size>::change_sign() const {
        return Physica::change_sign<Flags...>(toMachine());
    }

    template<Scalar T, int Size>
    SIMD<T, Size>& SIMD<T, Size>::cutoff(int count) {
        assert(0 < count && count < int(Size) && "[Error]: Invalid count");
        pack.cutoff(count);
        return *this;
    }

    template<Scalar T, int Size>
    T SIMD<T, Size>::sum() const {
        return Physica::horizontal_add(toMachine());
    }

    template<Scalar T, int Size>
    T SIMD<T, Size>::max() const {
        return Physica::horizontal_max1(toMachine());
    }

    template<Scalar T, int Size>
    T SIMD<T, Size>::min() const {
        return Physica::horizontal_min1(toMachine());
    }

    template<Scalar T, int Size>
    auto SIMD<T, Size>::isZero() const noexcept -> BoolSIMDType {
        return operator==(This(0));
    }

    template<Scalar T, int Size>
    auto SIMD<T, Size>::isSubNormal() const noexcept -> BoolSIMDType {
        return BoolSIMDType(Physica::is_zero_or_subnormal(pack));
    }

    template<Scalar T, int Size>
    auto SIMD<T, Size>::isPositive() const noexcept -> BoolSIMDType {
        return operator>(This(0));
    }

    template<Scalar T, int Size>
    auto SIMD<T, Size>::isNegative() const noexcept -> BoolSIMDType {
        return operator<(This(0));
    }

    template<Scalar T, int Size>
    auto SIMD<T, Size>::isFinite() const noexcept -> BoolSIMDType {
        return BoolSIMDType(is_finite(toMachine()));
    }

    template<Scalar T, int Size>
    auto SIMD<T, Size>::inf() noexcept -> SIMD {
        if constexpr (T::Prec == Float32) {
            if constexpr (Size == 4)
                return reinterpret_f(Vec4i(0x7F800000));
            else if constexpr (Size == 8)
                return reinterpret_f(Vec8i(0x7F800000));
            else {
                static_assert(Size == 16);
                return reinterpret_f(Vec16i(0x7F800000));
            }
        }
        else {
            if constexpr (Size == 2)
                return reinterpret_d(Vec2q(0x7FF0000000000000));
            else if constexpr (Size == 4)
                return reinterpret_d(Vec4q(0x7FF0000000000000));
            else {
                static_assert(Size == 8);
                return reinterpret_d(Vec8q(0x7FF0000000000000));
            }
        }
    }

    template<Scalar T, int Size>
    template<int... Order>
    SIMD<T, Size> SIMD<T, Size>::blend(const SIMD x, const SIMD y) {
        static_assert(sizeof...(Order) == Size, "[Error]: Size of Order do not match the packet");
        if constexpr (Size == 2)
            return Physica::blend2<Order...>(x.toMachine(), y.toMachine());
        else if constexpr (Size == 4)
            return Physica::blend4<Order...>(x.toMachine(), y.toMachine());
        else if constexpr (Size == 8)
            return Physica::blend8<Order...>(x.toMachine(), y.toMachine());
        else {
            static_assert(Size == 16, "[Error]: Unexpected size");
            return Physica::blend16<Order...>(x.toMachine(), y.toMachine());
        }
    }

    template<Scalar T, int Size>
    template<bool... Flags>
    auto SIMD<T, Size>::makeSignBits() noexcept -> This {
        using IntType = std::conditional<T::Prec == Float32, uint32_t, uint64_t>::type;
        constexpr auto Functor = [](bool flag) consteval static noexcept -> IntType {
            constexpr IntType Mask = T::Prec == Float32 ? IntType(0x80000000U) : IntType(0x8000000000000000U);
            return flag ? Mask : 0;
        };

        static_assert(sizeof...(Flags) == Size, "[Error]: Size do not match");
        if constexpr (T::Prec == Float32) {
            if constexpr (Size == 4)
                return __m128(_mm_setr_epi32(Functor(Flags)...));
            else if constexpr (Size == 8)
                return __m256(_mm256_setr_epi32(Functor(Flags)...));
            else {
                static_assert(Size == 16);
                // GCC 9.5.0 defines _mm512_setr_epi32 and _mm512_setr_epi64 as macro definitions, and we cannot apply template parameter pack expansion to macro definitions.
                // Therefore, it is necessary to define our own versions manually.
                const static auto mm512_setr_epi32 = [](int32_t i0, int32_t i1, int32_t i2, int32_t i3, int32_t i4, int32_t i5, int32_t i6, int32_t i7, int32_t i8, int32_t i9, int32_t i10, int32_t i11, int32_t i12, int32_t i13, int32_t i14, int32_t i15) static noexcept {
                    return _mm512_set_epi32(i15, i14, i13, i12, i11, i10, i9, i8, i7, i6, i5, i4, i3, i2, i1, i0);
                };
                return __m512(mm512_setr_epi32(Functor(Flags)...));
            }
        }
        else {
            if constexpr (Size == 2)
                return __m128d(_mm_setr_epi64(__m64(Functor(Flags))...));
            else if constexpr (Size == 4)
                return __m256d(_mm256_setr_epi64x(Functor(Flags)...));
            else {
                static_assert(Size == 8);
                const static auto mm512_setr_epi64 = [](int64_t i0, int64_t i1, int64_t i2, int64_t i3, int64_t i4, int64_t i5, int64_t i6, int64_t i7) static noexcept {
                    return _mm512_set_epi64(i7, i6, i5, i4, i3, i2, i1, i0);
                };
                return __m512d(mm512_setr_epi64(Functor(Flags)...));
            }
        }
    }

    template<Scalar T, int Size>
    template<RNG R>
    SIMD<T, Size> SIMD<T, Size>::random_uniform() {
        SIMD result{};
        std::array<T, Size> buffer{};
        for (auto& elem : buffer)
            elem = T::template random_uniform<R>();
        result.load(buffer.data());
        return result;
    }

    template<Scalar T, int Size>
    auto SIMD<T, Size>::zeros() noexcept -> SIMD {
        return SIMD(Pack::Zeros());
    }

    template<Scalar T, int Size>
    SIMD<T, Size> SIMD<T, Size>::select(BoolSIMDType flags, const SIMD x, const SIMD y) {
        return SIMD(Physica::select(flags.toMachine(), x.toMachine(), y.toMachine()));
    }

    template<Scalar T, int Size>
    template<int Order, int... Orders>
    constexpr unsigned int SIMD<T, Size>::makeShuffleMask(int order) {
        int result = order;
        if constexpr (T::Prec == Float32) {
            static_assert(0 <= Order && Order < 4, "[Error]: Invalid order");
            result |= Order << ((4 - (sizeof...(Orders) + 1)) * 2);
        }
        else {
            static_assert(0 <= Order && Order < 2, "[Error]: Invalid order");
            result |= Order << (Size - (sizeof...(Orders) + 1));
        }

        if constexpr (sizeof...(Orders) != 0)
            return makeShuffleMask<Orders...>(result);
        return result;
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<T, Size> nmul_add(const SIMD<T, Size> a, const SIMD<T, Size> b, const SIMD<T, Size> c) noexcept {
        static_assert(!T::isDiffable(), "[Error]: Not implemented");
        return SIMD<T, Size>(nmul_add(a.toMachine(), b.toMachine(), c.toMachine()));
    }

    template<Scalar T, int Size>
    [[nodiscard]] SIMD<T, Size> mul_sub(const SIMD<T, Size> a, const SIMD<T, Size> b, const SIMD<T, Size> c) noexcept {
        static_assert(!T::isDiffable(), "[Error]: Not implemented");
        return SIMD<T, Size>(mul_sub(a.toMachine(), b.toMachine(), c.toMachine()));
    }
    /**
     * \returns a * b -/+ c
     */
    template<Scalar T, int Size>
    [[nodiscard]] SIMD<T, Size> mul_addsub(const SIMD<T, Size> a, const SIMD<T, Size> b, const SIMD<T, Size> c) noexcept {
        static_assert(!T::isDiffable(), "[Error]: Not implemented");
        if constexpr (T::Prec == Float32) {
            if constexpr (Size == 4) {
                if constexpr (Instruset::hasFMA() || Instruset::hasAVX2())
                    return _mm_fmaddsub_ps(a, b, c);
                else if constexpr (Instruset::hasFMA4())
                    return _mm_maddsub_ps(a, b, c);
                else if constexpr (Instruset::hasSSE3())
                    return _mm_addsub_ps(a * b, c);
                else
                    return a * b + c.template change_sign<1, 0, 1, 0>();
            }
            else if constexpr (Size == 8) {
                if constexpr (Instruset::hasAVX()) {
                    if constexpr (Instruset::hasFMA() || Instruset::hasAVX2())
                        return _mm256_fmaddsub_ps(a, b, c);
                    else
                        return _mm256_addsub_ps(a * b, c);
                }
                else
                    return {mul_addsub(a.getLow(), b.getLow(), c.getLow()), mul_addsub(a.getHigh(), b.getHigh(), c.getHigh())};
            }
            else {
                static_assert(Size == 16, "[Error]: Unexpected size");
                if constexpr (Instruset::hasAVX512())
                    return _mm512_fmaddsub_ps(a, b, c);
                else
                    return {mul_addsub(a.getLow(), b.getLow(), c.getLow()), mul_addsub(a.getHigh(), b.getHigh(), c.getHigh())};
            }
        }
        else {
            static_assert(T::Prec == Float64, "[Error]: Unexpected float type");
            if constexpr (Size == 2) {
                if constexpr (Instruset::hasFMA() || Instruset::hasAVX2())
                    return _mm_fmaddsub_pd(a, b, c);
                else if constexpr (Instruset::hasFMA4())
                    return _mm_maddsub_pd(a, b, c);
                else if constexpr (Instruset::hasSSE3())
                    return _mm_addsub_pd(a * b, c);
                else
                    return a * b + c.template change_sign<1, 0>();
            }
            else if constexpr (Size == 4) {
                if constexpr (Instruset::hasAVX()) {
                    if constexpr (Instruset::hasFMA() || Instruset::hasAVX2())
                        return _mm256_fmaddsub_pd(a, b, c);
                    else
                        return _mm256_addsub_pd(a * b, c);
                }
                else
                    return {mul_addsub(a.getLow(), b.getLow(), c.getLow()), mul_addsub(a.getHigh(), b.getHigh(), c.getHigh())};
            }
            else {
                static_assert(Size == 8, "[Error]: Unexpected size");
                if constexpr (Instruset::hasAVX512())
                    return _mm512_fmaddsub_pd(a, b, c);
                else
                    return {mul_addsub(a.getLow(), b.getLow(), c.getLow()), mul_addsub(a.getHigh(), b.getHigh(), c.getHigh())};
            }
        }
    }
}

#include "Math.h"
