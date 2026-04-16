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

#include "../LValueVector.h"

namespace Physica {
    template<class Derived> class CompactVector;

    template<Vector V, size_t Length>
    class CompactVectorBlock : public CompactVector<CompactVectorBlock<V, Length>> {
        using This = CompactVectorBlock<V, Length>;
        using Base = CompactVector<This>;
    private:
        LazyDestroy<V> vec;
        size_t from;
        size_t to;
    public:
        CompactVectorBlock(V&& vec_, size_t from_, size_t to_);
        CompactVectorBlock(V&& vec_, size_t from_);
        CompactVectorBlock(const This& block) = default;
        CompactVectorBlock(This&&) noexcept = default;
        ~CompactVectorBlock() = default;
        /* Operators */
        using Base::operator=;
        This& operator=(const This& v);
        This& operator=(This&& v) noexcept;
        /* Operations */
        template<int Size> [[nodiscard]] auto packet(size_t index) const noexcept;
        template<int Size> [[nodiscard]] auto packet(size_t index, size_t count) const noexcept;
        void writePacket(Packet auto packet, size_t index) noexcept;
        void writePacket(Packet auto packet, size_t index, size_t count) noexcept;

        template<size_t Length_ = Dynamic>
        [[nodiscard]] auto head(this auto&&, size_t to = Length_) noexcept;
        template<size_t Length_ = Dynamic>
        [[nodiscard]] auto tail(this auto&&, size_t from) noexcept;
        template<size_t Length_ = Dynamic>
        [[nodiscard]] auto segment(this auto&&, size_t from, size_t to) noexcept;

        using Base::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }

        [[nodiscard]] auto values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] auto grads(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept;
        [[nodiscard]] auto data(this auto&&) noexcept;
    };

    template<Vector V, size_t Length>
    CompactVectorBlock<V, Length>::CompactVectorBlock(V&& vec_, size_t from_, size_t to_)
            : vec(std::forward<V>(vec_)), from(from_), to(to_) {
        assert(from_ < to);
        assert(to <= vec.getLength());
        assert(Length == Dynamic || Length == getLength());
    }

    template<Vector V, size_t Length>
    CompactVectorBlock<V, Length>::CompactVectorBlock(V&& vec, size_t from_)
            : CompactVectorBlock(std::forward<V>(vec), from_, vec.getLength()) {}

    template<Vector V, size_t Length>
    auto CompactVectorBlock<V, Length>::operator=(const This& v) -> This& {
        Base::operator=(v);
        return *this;
    }

    template<Vector V, size_t Length>
    auto CompactVectorBlock<V, Length>::operator=(This&& v) noexcept -> This& {
        Base::operator=(v);
        return *this;
    }

    template<Vector V, size_t Length>
    template<int Size>
    auto CompactVectorBlock<V, Length>::packet(size_t index) const noexcept {
        return vec.template packet<Size>(from + index);
    }

    template<Vector V, size_t Length>
    template<int Size>
    auto CompactVectorBlock<V, Length>::packet(size_t index, size_t count) const noexcept {
        return vec.template packet<Size>(from + index, count);
    }

    template<Vector V, size_t Length>
    void CompactVectorBlock<V, Length>::writePacket(Packet auto packet, size_t index) noexcept {
        return vec.writePacket(packet, from + index);
    }

    template<Vector V, size_t Length>
    void CompactVectorBlock<V, Length>::writePacket(Packet auto packet, size_t index, size_t count) noexcept {
        return vec.writePacket(packet, from + index, count);
    }

    template<Vector V, size_t Length>
    template<size_t Length_>
    auto CompactVectorBlock<V, Length>::head(this auto&& self, size_t to) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), V>(self.vec);
        using V1 = decltype(v);
        return CompactVectorBlock<V1, Length_>(std::forward<V1>(v), self.from, self.from + to);
    }

    template<Vector V, size_t Length>
    template<size_t Length_>
    auto CompactVectorBlock<V, Length>::tail(this auto&& self, size_t from) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), V>(self.vec);
        using V1 = decltype(v);
        return CompactVectorBlock<V1, Length_>(std::forward<V1>(v), self.from + from, self.to);
    }

    template<Vector V, size_t Length>
    template<size_t Length_>
    auto CompactVectorBlock<V, Length>::segment(this auto&& self, size_t from, size_t to) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), V>(self.vec);
        using V1 = decltype(v);
        return CompactVectorBlock<V1, Length_>(std::forward<V1>(v), self.from + from, self.from + to);
    }

    template<Vector V, size_t Length>
    auto CompactVectorBlock<V, Length>::values(this auto&& self) noexcept {
        auto&& v = propagate_rvalue_reference<decltype(self), V>(self.vec).values();
        using V1 = decltype(v);
        return CompactVectorBlock<V1, Length>(std::forward<V1>(v), self.from, self.to);
    }

    template<Vector V, size_t Length>
    template<int GradOrder>
    auto CompactVectorBlock<V, Length>::grads(this auto&& self) noexcept {
        auto&& g = propagate_rvalue_reference<decltype(self), V>(self.vec).template grads<GradOrder>();
        using V1 = decltype(g);
        return CompactVectorBlock<V1, Length>(std::forward<V1>(g), self.from, self.to);
    }

    template<Vector V, size_t Length>
    size_t CompactVectorBlock<V, Length>::getLength() const noexcept {
        if constexpr (Length == Dynamic)
            return to - from;
        else
            return Length;
    }

    template<Vector V, size_t Length>
    auto CompactVectorBlock<V, Length>::data(this auto&& self) noexcept {
        return self.vec.data() + self.from;
    }
}

namespace Physica {
    template<Vector V, size_t Length>
    class Traits<CompactVectorBlock<V, Length>> {
    public:
        using ScalarType = std::remove_cvref_t<V>::ScalarType;
        constexpr static size_t SizeAtCompile = Length;
        constexpr static bool FastAssign = false;

        using ElemType = ScalarType;
    };
}
