/*
 * Copyright 2021-2026 Weibo He.
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
        V& vec;
        size_t from;
        size_t to;
    public:
        LVectorBlock(V& vec_, size_t from_, size_t to_);
        LVectorBlock(V& vec_, size_t from_);
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
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept;
        [[nodiscard]] auto data_ptr(this auto&&, size_t index) noexcept;
    };

    template<Vector V, size_t Length>
    LVectorBlock<V, Length>::LVectorBlock(V& vec_, size_t from_, size_t to_)
            : vec(vec_), from(from_), to(to_) {
        assert(from_ < to);
        assert(to <= vec.getLength());
        assert(Length == Dynamic || Length == getLength());
    }

    template<Vector V, size_t Length>
    LVectorBlock<V, Length>::LVectorBlock(V& vec_, size_t from_) : LVectorBlock(vec_, from_, vec_.getLength()) {}

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
        return LVectorBlock<V, Length_>(self.vec, self.from, self.from + to);
    }

    template<Vector V, size_t Length>
    template<size_t Length_>
    auto LVectorBlock<V, Length>::tail(this auto&& self, size_t from) noexcept {
        return LVectorBlock<V, Length_>(self.vec, self.from + from, self.to);
    }

    template<Vector V, size_t Length>
    template<size_t Length_>
    auto LVectorBlock<V, Length>::segment(this auto&& self, size_t from, size_t to) noexcept {
        return LVectorBlock<V, Length_>(self.vec, self.from + from, self.from + to);
    }

    template<Vector V, size_t Length>
    size_t LVectorBlock<V, Length>::getLength() const noexcept {
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
}

namespace Physica {
    template<Vector V, size_t Length>
    class Traits<LVectorBlock<V, Length>> {
    public:
        using ScalarType = V::ScalarType;
        constexpr static size_t SizeAtCompile = Length;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };
}
