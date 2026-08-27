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
#include <cassert>
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Scalar/Real.h"
#include "Physica/Core/Scalar/RealImpl/SIMD.h"

using namespace Physica;
/**
 * Adapted from [1], we require a more strict range of \param n. SIMD<T, Size>::store() is similar.
 *
 * Partial store and partial load are typically used for tail elements, and it is less interesting to inline them.
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

#if INSTRSET >= 2
template void SIMD<float32, 4>::load(const float32*, int) & noexcept;
template void SIMD<float64, 2>::load(const float64*, int) & noexcept;
template void SIMD<float32, 4>::store(float32*, int) const noexcept;
template void SIMD<float64, 2>::store(float64*, int) const noexcept;
#endif

#if INSTRSET >= 7
template void SIMD<float32, 8>::load(const float32*, int) & noexcept;
template void SIMD<float64, 4>::load(const float64*, int) & noexcept;
template void SIMD<float32, 8>::store(float32*, int) const noexcept;
template void SIMD<float64, 4>::store(float64*, int) const noexcept;
#endif

#if INSTRSET >= 9
template void SIMD<float32, 16>::load(const float32*, int) & noexcept;
template void SIMD<float64, 8>::load(const float64*, int) & noexcept;
template void SIMD<float32, 16>::store(float32*, int) const noexcept;
template void SIMD<float64, 8>::store(float64*, int) const noexcept;
#endif
