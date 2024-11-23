/*
 * Copyright 2023-2024 Weibo He.
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

#include "Physica/Core/MultiPrecision/ExprType.h"
#include "Physica/Core/Utils/Unreachable.h"

namespace Physica::Core {
    template<Scalar T, size_t Size>
    SIMD<T, Size>::SIMD(HalfType a, HalfType b) : Base(a.toMachine(), b.toMachine()) {}

    template<Scalar T, size_t Size>
    SIMD<T, Size>::SIMD(T s, int count) {
        assert(0 < count && count <= int(Size) && "[Error]: Invalid count");
        if (count == Size) {
            *this = SIMD(s);
            return;
        }

        constexpr int HalfSize = Size / 2;
        if constexpr (isSeparatable) {
            if (count > HalfSize)
                *this = SIMD(HalfType(s), HalfType(s, count - HalfSize));
            else
                *this = SIMD(HalfType(s, count), HalfType(0));
        }
        else {
            if constexpr (T::Option == Float32) {
                const float f = float(s);
                switch(count) {
                case 1:
                    *this = SIMD(Base(f, 0.0F, 0.0F, 0.0F));
                    break;
                case 2:
                    *this = SIMD(Base(f, f, 0.0F, 0.0F));
                    break;
                case 3:
                    *this = SIMD(Base(f, f, f, 0.0F));
                    break;
                default:
                    unreachable();
                }
            }
            else {
                static_assert(T::Option == Float64, "[Error]: Unsupported float type");
                *this = SIMD(Base(double(s), 0.0));
            }
        }
    }

    template<Scalar T, size_t Size>
    inline T SIMD<T, Size>::operator[](int index) const {
        if constexpr (isForward)
            return T(Base::operator[](index * 2), Base::operator[](index * 2 + 1));
        else
            return T(Base::operator[](index));
    }

    template<Scalar T, size_t Size>
    inline SIMD<T, Size> SIMD<T, Size>::operator+(const SIMD& other) const {
        return SIMD(toMachine() + other.toMachine());
    }

    template<Scalar T, size_t Size>
    inline SIMD<T, Size> SIMD<T, Size>::operator-(const SIMD& other) const {
        return SIMD(toMachine() - other.toMachine());
    }

    template<Scalar T, size_t Size>
    inline SIMD<T, Size> SIMD<T, Size>::operator*(const SIMD& other) const {
        return SIMD(toMachine() * other.toMachine());
    }

    template<Scalar T, size_t Size>
    inline SIMD<T, Size> SIMD<T, Size>::operator*(const T& x) const {
        return operator*(SIMD(x));
    }

    template<Scalar T, size_t Size>
    inline SIMD<T, Size> SIMD<T, Size>::operator/(const SIMD& other) const {
        return SIMD(toMachine() / other.toMachine());
    }

    template<Scalar T, size_t Size>
    inline SIMD<T, Size> SIMD<T, Size>::operator-() const {
        return SIMD(-toMachine());
    }

    template<Scalar T, size_t Size>
    inline SIMD<T, Size> SIMD<T, Size>::operator&(const SIMD& other) const {
        return SIMD(toMachine() & other.toMachine());
    }

    template<Scalar T, size_t Size>
    inline SIMD<T, Size> SIMD<T, Size>::operator|(const SIMD& other) const {
        return SIMD(toMachine() | other.toMachine());
    }

    template<Scalar T, size_t Size>
    inline SIMD<T, Size> SIMD<T, Size>::operator^(const SIMD& other) const {
        return SIMD(toMachine() ^ other.toMachine());
    }

    template<Scalar T, size_t Size>
    inline SIMD<T, Size> SIMD<T, Size>::operator!() const {
        return SIMD(!toMachine());
    }

    template<Scalar T, size_t Size>
    inline auto SIMD<T, Size>::operator==(const SIMD& other) const {
        return BoolSIMDType(toMachine() == other.toMachine());
    }

    template<Scalar T, size_t Size>
    inline auto SIMD<T, Size>::operator>(const SIMD& other) const {
        return BoolSIMDType(toMachine() > other.toMachine());
    }
    
    template<Scalar T, size_t Size>
    inline auto SIMD<T, Size>::operator<(const SIMD& other) const {
        return BoolSIMDType(toMachine() < other.toMachine());
    }

    template<Scalar T, size_t Size>
    inline void SIMD<T, Size>::load(const T* p) {
        Base::load(reinterpret_cast<const MachineType*>(p));
    }

    template<Scalar T, size_t Size>
    inline void SIMD<T, Size>::load_partial(const T* p, int n) {
        Base::load_partial(n, reinterpret_cast<const MachineType*>(p));
    }

    template<Scalar T, size_t Size>
    inline void SIMD<T, Size>::store(T* p) const {
        Base::store(reinterpret_cast<MachineType*>(p));
    }

    template<Scalar T, size_t Size>
    inline void SIMD<T, Size>::store_partial(T* p, int n) const {
        Base::store_partial(n, reinterpret_cast<MachineType*>(p));
    }

    template<Scalar T, size_t Size>
    inline void SIMD<T, Size>::insert(int index, const T& value) {
        if constexpr (isForward) {
            Base::insert(index * 2, value.getValue().toMachine());
            Base::insert(index * 2 + 1, value.getGrad().toMachine());
        }
        else
            Base::insert(index, value.toMachine());
    }

    template<Scalar T, size_t Size>
    template<int... Order>
    inline SIMD<T, Size> SIMD<T, Size>::shuffle() const {
        constexpr unsigned int mask = makeShuffleMask<Order...>(0);
        if constexpr (T::Option == Float32) {
            static_assert(sizeof...(Order) == 4, "[Error]: Invalid number of orders");
            if constexpr (Size == 4)
                return _mm_shuffle_ps(*this, *this, mask);
            else if constexpr (Size == 8)
                return _mm256_shuffle_ps(*this, *this, mask);
            else if constexpr (Size == 16)
                return _mm512_shuffle_ps(*this, *this, mask);
        }
        else {
            static_assert(sizeof...(Order) == Size, "[Error]: Invalid number of orders");
            if constexpr (Size == 2)
                return _mm_shuffle_pd(*this, *this, mask);
            else if constexpr (Size == 4)
                return _mm256_shuffle_pd(*this, *this, mask);
            else if constexpr (Size == 8)
                return _mm512_shuffle_pd(*this, *this, mask);
        }
    }

    template<Scalar T, size_t Size>
    template<int... Order>
    inline SIMD<T, Size> SIMD<T, Size>::permute() const {
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

    template<Scalar T, size_t Size>
    template<int... Flags>
    inline SIMD<T, Size> SIMD<T, Size>::change_sign() const {
        return Physica::change_sign<Flags...>(toMachine());
    }

    template<Scalar T, size_t Size>
    inline SIMD<T, Size>& SIMD<T, Size>::cutoff(int count) {
        assert(0 < count && count < int(Size) && "[Error]: Invalid count");
        Base::cutoff(count);
        return *this;
    }

    template<Scalar T, size_t Size>
    inline T SIMD<T, Size>::sum() const {
        return Physica::horizontal_add(toMachine());
    }

    template<Scalar T, size_t Size>
    inline T SIMD<T, Size>::max() const {
        return Physica::horizontal_max(toMachine());
    }

    template<Scalar T, size_t Size>
    inline T SIMD<T, Size>::min() const {
        return Physica::horizontal_min(toMachine());
    }

    template<Scalar T, size_t Size>
    template<bool... Flags>
    SIMD<T, Size> SIMD<T, Size>::makeSignBits() {
        using IntType = std::conditional<T::Option == Float32, uint32_t, uint64_t>::type;
        constexpr auto Functor = [](bool flag) constexpr -> IntType {
            constexpr IntType Mask = T::Option == Float32 ? IntType(0x80000000U) : IntType(0x8000000000000000U);
            return flag ? Mask : 0;
        };

        static_assert(sizeof...(Flags) == Size, "[Error]: Size do not match");
        if constexpr (T::Option == Float32) {
            if constexpr (Size == 4)
                return __m128(_mm_setr_epi32(Functor(Flags)...));
            else if constexpr (Size == 8)
                return __m256(_mm256_setr_epi32(Functor(Flags)...));
            else {
                static_assert(Size == 16);
                // GCC 9.5.0 defines _mm512_setr_epi32 and _mm512_setr_epi64 as macro definitions, and we cannot apply template parameter pack expansion to macro definitions.
                // Therefore, it is necessary to define our own versions manually.
                const static auto mm512_setr_epi32 = [](int32_t i0, int32_t i1, int32_t i2, int32_t i3, int32_t i4, int32_t i5, int32_t i6, int32_t i7,
                                                        int32_t i8, int32_t i9, int32_t i10, int32_t i11, int32_t i12, int32_t i13, int32_t i14, int32_t i15) {
                    return _mm512_set_epi32(i15, i14, i13, i12, i11, i10, i9, i8, i7, i6, i5, i4, i3, i2, i1, i0);
                };
                return __m512(mm512_setr_epi32(Functor(Flags)...));
            }
        }
        else {
            if constexpr (Size == 2)
                return __m128d(_mm_setr_epi64(Functor(Flags)...));
            else if constexpr (Size == 4)
                return __m256d(_mm256_setr_epi64(Functor(Flags)...));
            else {
                static_assert(Size == 8);
                const static auto mm512_setr_epi64 = [](int64_t i0, int64_t i1, int64_t i2, int64_t i3, int64_t i4, int64_t i5, int64_t i6, int64_t i7) {
                    return _mm512_set_epi64(i7, i6, i5, i4, i3, i2, i1, i0);
                };
                return __m512d(mm512_setr_epi64(Functor(Flags)...));
            }
        }
    }

    template<Scalar T, size_t Size>
    template<RandomGenerator R>
    SIMD<T, Size> SIMD<T, Size>::random_uniform() {
        SIMD result{};
        T buffer[Size];
        for (auto& elem : buffer)
            elem = T::template random_uniform<R>();
        result.load(buffer);
        return result;
    }

    template<Scalar T, size_t Size>
    SIMD<T, Size> SIMD<T, Size>::select(BoolSIMDType flags, const SIMD& x, const SIMD& y) {
        return SIMD(Physica::select(flags.toMachine(), x.toMachine(), y.toMachine()));
    }

    template<Scalar T, size_t Size>
    template<int Order, int... Orders>
    constexpr unsigned int SIMD<T, Size>::makeShuffleMask(int order) {
        int result = order;
        if constexpr (T::Option == Float32) {
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
    //////////////////////////////////////////////////////////////////
    template<Scalar T, size_t Size>
    [[nodiscard]] inline SIMD<T, Size> mul_add(
            const SIMD<T, Size>& a,
            const SIMD<T, Size>& b,
            const SIMD<T, Size>& c) {
        return SIMD<T, Size>(mul_add(a.toMachine(), b.toMachine(), c.toMachine()));
    }

    template<Scalar T, size_t Size>
    [[nodiscard]] inline SIMD<T, Size> nmul_add(
            const SIMD<T, Size> a,
            const SIMD<T, Size> b,
            const SIMD<T, Size> c) {
        static_assert(!T::isDifferentiable, "[Error]: Not implemented");
        return SIMD<T, Size>(nmul_add(a.toMachine(), b.toMachine(), c.toMachine()));
    }

    template<Scalar T, size_t Size>
    [[nodiscard]] inline SIMD<T, Size> mul_sub(
            const SIMD<T, Size> a,
            const SIMD<T, Size> b,
            const SIMD<T, Size> c) {
        static_assert(!T::isDifferentiable, "[Error]: Not implemented");
        return SIMD<T, Size>(mul_sub(a.toMachine(), b.toMachine(), c.toMachine()));
    }

    template<Scalar T, size_t Size>
    [[nodiscard]] inline SIMD<T, Size> mul_addsub(
            const SIMD<T, Size> a,
            const SIMD<T, Size> b,
            const SIMD<T, Size> c) {
        static_assert(!T::isDifferentiable, "[Error]: Not implemented");
        if constexpr (T::Option == Float32) {
            if constexpr (Size == 4) {
                if constexpr (Instruset::hasFMA() || Instruset::hasAVX2())
                    return  _mm_fmaddsub_ps(a, b, c);
                else if constexpr (Instruset::hasFMA4())
                    return  _mm_maddsub_ps(a, b, c);
                else if constexpr (Instruset::hasSSE3())
                    return _mm_addsub_ps(a * b, c);
                else
                    return a * b + c.template change_sign<1, 0, 1, 0>();
            }
            else if constexpr (Size == 8) {
                if constexpr (Instruset::hasAVX()) {
                    if constexpr (Instruset::hasFMA() || Instruset::hasAVX2())
                        return  _mm256_fmaddsub_ps(a, b, c);
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
            static_assert(T::Option == Float64, "[Error]: Unexpected float type");
            if constexpr (Size == 2) {
                if constexpr (Instruset::hasFMA() || Instruset::hasAVX2())
                    return  _mm_fmaddsub_pd(a, b, c);
                else if constexpr (Instruset::hasFMA4())
                    return  _mm_maddsub_pd(a, b, c);
                else if constexpr (Instruset::hasSSE3())
                    return _mm_addsub_pd(a * b, c);
                else
                    return a * b + c.template change_sign<1, 0>();
            }
            else if constexpr (Size == 4) {
                if constexpr (Instruset::hasAVX()) {
                    if constexpr (Instruset::hasFMA() || Instruset::hasAVX2())
                        return  _mm256_fmaddsub_pd(a, b, c);
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
