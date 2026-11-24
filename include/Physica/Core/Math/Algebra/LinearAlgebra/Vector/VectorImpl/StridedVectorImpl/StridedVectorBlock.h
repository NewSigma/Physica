/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "../LValueVector.h"

namespace Physica {
    template<class Derived> class StridedVector;

    template<Vector V, size_t Length>
    class StridedVectorBlock : public StridedVector<StridedVectorBlock<V, Length>> {
        using This = StridedVectorBlock<V, Length>;
        using Base = StridedVector<This>;
    private:
        decay_rvalue_t<V> vec;
        size_t from;
        size_t to;
    public:
        constexpr StridedVectorBlock(V&& vec_, size_t from, size_t to);
        constexpr StridedVectorBlock(V&& vec_, size_t from);
        StridedVectorBlock(const This& block) = default;
        StridedVectorBlock(This&&) noexcept = default;
        ~StridedVectorBlock() = default;
        /* Operators */
        using Base::operator=;
        This& operator=(const This& v);
        This& operator=(This&& v) noexcept;
        /* Operations */
        using Base::resize;
        constexpr void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] constexpr size_t getLength() const noexcept;
        [[nodiscard]] constexpr size_t getStride() const noexcept;
        [[nodiscard]] constexpr auto data_handle(this auto&& self) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getStrideAtCompile() noexcept;
    };

    template<Vector V, size_t Length>
    constexpr StridedVectorBlock<V, Length>::StridedVectorBlock(V&& vec_, size_t from, size_t to)
            : vec(std::forward<V>(vec_)), from(from), to(to) {
        assert(from < to);
        assert(to <= vec.getLength());
        assert(Length == Dynamic || Length == getLength());
    }

    template<Vector V, size_t Length>
    constexpr StridedVectorBlock<V, Length>::StridedVectorBlock(V&& vec_, size_t from_)
            : StridedVectorBlock(std::forward<V>(vec_), from_, vec_.getLength()) {}

    template<Vector V, size_t Length>
    auto StridedVectorBlock<V, Length>::operator=(const This& v) -> This& {
        Base::operator=(v);
        return *this;
    }

    template<Vector V, size_t Length>
    auto StridedVectorBlock<V, Length>::operator=(This&& v) noexcept -> This& {
        Base::operator=(v);
        return *this;
    }

    template<Vector V, size_t Length>
    auto StridedVectorBlock<V, Length>::values(this auto&& self) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), V>(self.vec).values();
        using V1 = decltype(v);
        if constexpr (std::remove_cvref_t<V1>::isStrided())
            return StridedVectorBlock<V1, Length>(std::forward<V1>(v), self.from, self.to);
        else if constexpr (std::remove_cvref_t<V1>::isLValueVector())
            return LVectorBlock<V1, Length>(std::forward<V1>(v), self.from, self.to);
        else
            return RVectorBlock<V1, Length>(std::forward<V1>(v), self.from, self.to);
    }

    template<Vector V, size_t Length>
    constexpr size_t StridedVectorBlock<V, Length>::getLength() const noexcept {
        if constexpr (Length == Dynamic)
            return to - from;
        else
            return Length;
    }

    template<Vector V, size_t Length>
    constexpr size_t StridedVectorBlock<V, Length>::getStride() const noexcept {
        return vec.getStride();
    }

    template<Vector V, size_t Length>
    constexpr auto StridedVectorBlock<V, Length>::data_handle(this auto&& self) noexcept {
        return self.vec.data_handle() + self.from * self.vec.getStride();
    }

    template<Vector V, size_t Length>
    __host__ __device__ consteval size_t StridedVectorBlock<V, Length>::getSizeAtCompile() noexcept {
        return Length;
    }

    template<Vector V, size_t Length>
    __host__ __device__ consteval size_t StridedVectorBlock<V, Length>::getStrideAtCompile() noexcept {
        return std::remove_cvref_t<V>::getStrideAtCompile();
    }
}

namespace Physica {
    template<Vector V, size_t Length>
    class Traits<StridedVectorBlock<V, Length>> {
    public:
        using ScalarType = std::remove_cvref_t<V>::ScalarType;
    };
}
