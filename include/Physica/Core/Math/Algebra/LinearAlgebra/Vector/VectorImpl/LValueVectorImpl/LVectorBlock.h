/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "../LValueVector.h"

namespace Physica {
    template<class Derived> class LValueVector;
    /**
     * Reference a part of the given vector
     */
    template<Vector V, size_t Length>
    class LVectorBlock : public LValueVector<LVectorBlock<V, Length>> {
        using This = LVectorBlock<V, Length>;
        using Base = LValueVector<This>;
    public:
        using ScalarType = Base::ScalarType;
    private:
        decay_rvalue_t<V> vec;
        size_t from;
        size_t to;
    public:
        LVectorBlock(V&& vec_, size_t from, size_t to);
        LVectorBlock(V&& vec_, size_t from);
        LVectorBlock(const This& block) = default;
        LVectorBlock(This&&) noexcept = default;
        ~LVectorBlock() = default;
        /* Operators */
        using Base::operator=;
        This& operator=(const This& v);
        This& operator=(This&& v) noexcept;
        /* Operations */
        template<size_t Length_ = Dynamic>
        [[nodiscard]] auto head(this auto&&, size_t to = Length_) noexcept;
        template<size_t Length_ = Dynamic>
        [[nodiscard]] auto tail(this auto&&, size_t from) noexcept;
        template<size_t Length_ = Dynamic>
        [[nodiscard]] auto segment(this auto&&, size_t from, size_t to) noexcept;

        using Base::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] constexpr size_t getLength() const noexcept;
        [[nodiscard]] auto data_ptr(this auto&&, size_t index) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept;
    };

    template<Vector V, size_t Length>
    LVectorBlock<V, Length>::LVectorBlock(V&& vec_, size_t from, size_t to)
            : vec(std::forward<V>(vec_)), from(from), to(to) {
        assert(from < to);
        assert(to <= vec.getLength());
        assert(Length == Dynamic || Length == getLength());
    }

    template<Vector V, size_t Length>
    LVectorBlock<V, Length>::LVectorBlock(V&& vec_, size_t from_) : LVectorBlock(std::forward<V>(vec_), from_, vec_.getLength()) {}

    template<Vector V, size_t Length>
    auto LVectorBlock<V, Length>::operator=(const This& v) -> This& {
        Base::operator=(v);
        return *this;
    }

    template<Vector V, size_t Length>
    auto LVectorBlock<V, Length>::operator=(This&& v) noexcept -> This& {
        Base::operator=(v);
        return *this;
    }

    template<Vector V, size_t Length>
    template<size_t Length_>
    auto LVectorBlock<V, Length>::head(this auto&& self, size_t to) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), V>(self.vec);
        using V1 = decltype(v);
        return LVectorBlock<V1, Length_>(std::forward<V1>(v), self.from, self.from + to);
    }

    template<Vector V, size_t Length>
    template<size_t Length_>
    auto LVectorBlock<V, Length>::tail(this auto&& self, size_t from) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), V>(self.vec);
        using V1 = decltype(v);
        return LVectorBlock<V1, Length_>(std::forward<V1>(v), self.from + from, self.to);
    }

    template<Vector V, size_t Length>
    template<size_t Length_>
    auto LVectorBlock<V, Length>::segment(this auto&& self, size_t from, size_t to) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), V>(self.vec);
        using V1 = decltype(v);
        return LVectorBlock<V1, Length_>(std::forward<V1>(v), self.from + from, self.from + to);
    }

    template<Vector V, size_t Length>
    auto LVectorBlock<V, Length>::values(this auto&& self) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), V>(self.vec).values();
        using V1 = decltype(v);
        if constexpr (v.isLValueVector())
            return LVectorBlock<V1, Length>(std::forward<V1>(v), self.from, self.to);
        else
            return RVectorBlock<V1, Length>(std::forward<V1>(v), self.from, self.to);
    }

    template<Vector V, size_t Length>
    constexpr size_t LVectorBlock<V, Length>::getLength() const noexcept {
        if constexpr (Length == Dynamic)
            return to - from;
        else
            return Length;
    }

    template<Vector V, size_t Length>
    auto LVectorBlock<V, Length>::data_ptr(this auto&& self, size_t index) noexcept {
        assert((self.from + index) < self.to);
        return self.vec.data_ptr(index + self.from);
    }

    template<Vector V, size_t Length>
    __host__ __device__ consteval size_t LVectorBlock<V, Length>::getSizeAtCompile() noexcept {
        return Length;
    }
}

namespace Physica {
    template<Vector V, size_t Length>
    class Traits<LVectorBlock<V, Length>> {
        static_assert(std::remove_cvref_t<V>::isLValueVector());
    public:
        using ScalarType = std::remove_cvref_t<V>::ScalarType;
    };
}
