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
    template<class Derived> class ContinuousVector;

    template<Vector V, size_t Length>
    class ContinuousVectorBlock : public ContinuousVector<ContinuousVectorBlock<V, Length>> {
        using This = ContinuousVectorBlock<V, Length>;
        using Base = ContinuousVector<This>;
    private:
        V& vec;
        size_t from;
        size_t to;
    public:
        ContinuousVectorBlock(ContinuousVector<V>& vec_, size_t from_, size_t to_);
        ContinuousVectorBlock(ContinuousVector<V>& vec_, size_t from_);
        ContinuousVectorBlock(const This& block) = default;
        ContinuousVectorBlock(This&&) noexcept = default;
        ~ContinuousVectorBlock() = default;
        /* Operators */
        using Base::operator=;
        This& operator=(const This& v);
        This& operator=(This&& v) noexcept;
        /* Operations */
        template<Packet Pack> [[nodiscard]] Pack packet(size_t index) const;
        template<Packet Pack> [[nodiscard]] Pack packetPartial(size_t index, size_t count) const;
        void writePacket(size_t index, Packet auto packet);
        void writePacketPartial(size_t index, size_t count, Packet auto packet);

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
    ContinuousVectorBlock<V, Length>::ContinuousVectorBlock(ContinuousVector<V>& vec_, size_t from_, size_t to_)
            : vec(vec_.getDerived()), from(from_), to(to_) {
        assert(from_ < to);
        assert(to <= vec.getLength());
        assert(Length == Dynamic || Length == getLength());
    }

    template<Vector V, size_t Length>
    ContinuousVectorBlock<V, Length>::ContinuousVectorBlock(ContinuousVector<V>& vec_, size_t from_)
            : ContinuousVectorBlock(vec_, from_, vec_.getLength()) {}

    template<Vector V, size_t Length>
    auto ContinuousVectorBlock<V, Length>::operator=(const This& v) -> This& {
        Base::operator=(v);
        return *this;
    }

    template<Vector V, size_t Length>
    auto ContinuousVectorBlock<V, Length>::operator=(This&& v) noexcept -> This& {
        Base::operator=(v);
        return *this;
    }

    template<Vector V, size_t Length>
    template<Packet Pack>
    Pack ContinuousVectorBlock<V, Length>::packet(size_t index) const {
        return vec.template packet<Pack>(from + index);
    }

    template<Vector V, size_t Length>
    template<Packet Pack>
    Pack ContinuousVectorBlock<V, Length>::packetPartial(size_t index, size_t count) const {
        return vec.template packetPartial<Pack>(from + index, count);
    }

    template<Vector V, size_t Length>
    void ContinuousVectorBlock<V, Length>::writePacket(size_t index, const Packet auto packet) {
        return vec.writePacket(from + index, packet);
    }

    template<Vector V, size_t Length>
    void ContinuousVectorBlock<V, Length>::writePacketPartial(size_t index, size_t count, const Packet auto packet) {
        return vec.writePacketPartial(from + index, count, packet);
    }

    template<Vector V, size_t Length>
    auto ContinuousVectorBlock<V, Length>::values(this auto&& self) noexcept {
        auto&& v = self.vec.values();
        using V1 = std::remove_reference<decltype(v)>::type;
        return ContinuousVectorBlock<V1, Length>(v, self.from, self.to);
    }

    template<Vector V, size_t Length>
    template<int GradOrder>
    auto ContinuousVectorBlock<V, Length>::grads(this auto&& self) noexcept {
        auto&& g = self.vec.template grads<GradOrder>();
        using V1 = std::remove_reference<decltype(g)>::type;
        return ContinuousVectorBlock<V1, Length>(g, self.from, self.to);
    }

    template<Vector V, size_t Length>
    size_t ContinuousVectorBlock<V, Length>::getLength() const noexcept {
        if constexpr (Length == Dynamic)
            return to - from;
        else
            return Length;
    }

    template<Vector V, size_t Length>
    auto ContinuousVectorBlock<V, Length>::data(this auto&& self) noexcept {
        return self.vec.data() + self.from;
    }
}

namespace Physica {
    template<Vector V, size_t Length>
    class Traits<ContinuousVectorBlock<V, Length>> {
    public:
        using ScalarType = V::ScalarType;
        constexpr static size_t SizeAtCompile = Length;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = true;

        using ElemType = ScalarType;
    };
}
