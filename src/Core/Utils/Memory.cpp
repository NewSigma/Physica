/*
 * Copyright 2026 Weibo He.
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
#include <cstring>
#include <memory>
#include <type_traits>
#include "Physica/Core/Scalar/Real.h"
#include "Physica/Core/Utils/Unroller.h"
#include "Physica/Core/Utils/Builtin.h"
#include "Physica/Core/Utils/Memory.h"

#if defined(__x86_64__) || defined(_M_X64) || defined(__SSE2__)
    #include <immintrin.h>
    #include "Physica/Core/Scalar/RealImpl/SIMDImpl/Instruset.h"
    #define PHYSICA_HAS_STREAM_NT
#endif

using namespace Physica;

namespace {
    // Cached and non-temporal stores must not share a cache line
    constexpr size_t LineSizeL1D = HostDevAttr::LineSizeL1D;

    void memcpy_nt_nofence(void* dst, const void* src, size_t n) noexcept {
    #ifdef PHYSICA_HAS_STREAM_NT
        using Pack = std::conditional_t<Instruset::hasAVX512(), SIMD<float64, 8>,
                                        std::conditional_t<Instruset::hasAVX2(), SIMD<float32, 8>, SIMD<float32, 4>>>;
        using T = Pack::ScalarType;
        constexpr size_t ElemBytes = sizeof(T);
        void* aligned = dst;
        size_t space = n;
        if (std::align(LineSizeL1D, LineSizeL1D, aligned, space) == nullptr) {
            std::memcpy(dst, src, n);
            return;
        }

        auto* dst_bytes = static_cast<std::byte*>(dst);
        const auto* src_bytes = static_cast<const std::byte*>(src);
        auto* d = Physica::assume_aligned<LineSizeL1D>(static_cast<std::byte*>(aligned));
        const size_t head = d - dst_bytes;
        if (head)
            std::memcpy(dst_bytes, src_bytes, head);

        const size_t size_aligned = (n - head) & ~(LineSizeL1D - 1);
        const auto* s_aligned = src_bytes + head;
        Unroller<Pack, 2> unroller{};
        [[maybe_unused]] const size_t consumed = unroller.loop_recursive([d, s_aligned](Pack buffer, size_t index) noexcept {
            const size_t offset = index * ElemBytes;
            buffer.load(reinterpret_cast<const T*>(s_aligned + offset));
            buffer.store_nt(reinterpret_cast<T*>(d + offset));
            return buffer;
        }, size_aligned / ElemBytes);
        assert(d + consumed * ElemBytes == d + size_aligned && "[Error]: Unexpected size mismatch");

        if (size_aligned < n - head)
            std::memcpy(d + size_aligned, s_aligned + size_aligned, n - head - size_aligned);
    #else
        std::memcpy(dst, src, n);
    #endif
    }

    [[gnu::always_inline]] inline void sfence() noexcept {
    #ifdef PHYSICA_HAS_STREAM_NT
        _mm_sfence();
    #endif
    }
}

namespace Physica {
    void memcpy2D(void* dst, size_t dstride, const void* src, size_t sstride, size_t length, size_t numRow) noexcept {
        assert((dstride >= length) && (sstride >= length) && "[Error]: Storage overlaps");
        if (length == 0 || numRow == 0) [[unlikely]]
            return;

        const bool isCompact = dstride == length && sstride == length;
        if (isCompact) [[unlikely]] { // Prefer std::memcpy in caller
            std::memcpy(dst, src, length * numRow);
            return;
        }

        auto* d = static_cast<std::byte*>(dst);
        const auto* s = static_cast<const std::byte*>(src);
    #ifdef PHYSICA_HAS_STREAM_NT
        if (numRow * length >= HostDevAttr::ThresholdStreamNT && length >= LineSizeL1D) {
            for (size_t row = 0; row < numRow; ++row)
                memcpy_nt_nofence(d + row * dstride, s + row * sstride, length);
            sfence();
            return;
        }
    #endif

        if (length > LineSizeL1D) {
            for (size_t row = 0; row < numRow; ++row)
                std::memcpy(d + row * dstride, s + row * sstride, length);
        }
        else {
            for (size_t row = 0; row < numRow; ++row) {
                for (size_t i = 0; i < length; ++i)
                    d[row * dstride + i] = s[row * sstride + i];
            }
        }
    }

    void memcpy_nt(void* dst, const void* src, size_t n) noexcept {
        memcpy_nt_nofence(dst, src, n);
        sfence();
    }
}
