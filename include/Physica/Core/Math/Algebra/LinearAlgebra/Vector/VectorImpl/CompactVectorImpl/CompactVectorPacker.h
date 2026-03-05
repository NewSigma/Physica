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
#pragma once

#include "../CompactVector.h"

namespace Physica {
    template<Vector V, Packet Pack>
    class CompactVectorPacker : public LVectorPacker<V, Pack> {
        using This = CompactVectorPacker<V, Pack>;
        using Base = LVectorPacker<V, Pack>;
        using PtrTy = V::ScalarType::ConstPtrTy;
    private:
        PtrTy data;
    public:
        CompactVectorPacker(const V& vec);
        CompactVectorPacker(const This&) = default;
        CompactVectorPacker(This&&) noexcept = default;
        ~CompactVectorPacker() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        [[nodiscard]] Pack load(size_t index) const noexcept;
        [[nodiscard]] Pack load(size_t index, size_t count) const noexcept;
        void store(Pack pack, size_t index) noexcept;
        void store(Pack pack, size_t index, size_t count) noexcept;
    };

    template<Vector V, Packet Pack>
    CompactVectorPacker<V, Pack>::CompactVectorPacker(const V& vec) : Base(vec), data(vec.data()) {}

    template<Vector V, Packet Pack>
    auto CompactVectorPacker<V, Pack>::load(size_t index) const noexcept -> Pack {
        assert(index + Pack::size() <= Base::vector().getLength());
        Pack pack{};
        pack.load(data + index);
        return pack;
    }

    template<Vector V, Packet Pack>
    auto CompactVectorPacker<V, Pack>::load(size_t index, size_t count) const noexcept -> Pack {
        assert(index + count <= Base::vector().getLength());
        assert(0 < count && count < Pack::size() && "[Error]: Invalid size for partial operation");
        Pack pack{};
        pack.load(data + index, count);
        return pack;
    }

    template<Vector V, Packet Pack>
    void CompactVectorPacker<V, Pack>::store(const Pack packet, size_t index) noexcept {
        assert(index + Pack::size() <= Base::vector().getLength());
        packet.store(data + index);
    }

    template<Vector V, Packet Pack>
    void CompactVectorPacker<V, Pack>::store(const Pack packet, size_t index, size_t count) noexcept {
        assert(index + count <= Base::vector().getLength());
        assert(0 < count && count < Pack::size() && "[Error]: Invalid size for partial operation");
        packet.store(data + index, count);
    }
}
