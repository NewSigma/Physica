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

#include "../LValueVector.h"

namespace Physica {
    template<Vector V, Packet Pack>
    class LVectorPacker : public RVectorPacker<V, Pack> {
        using This = LVectorPacker<V, Pack>;
        using Base = RVectorPacker<V, Pack>;
    public:
        LVectorPacker(const V& vec) : Base(vec) {}
        LVectorPacker(const This&) = default;
        LVectorPacker(This&&) noexcept = default;
        ~LVectorPacker() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        void store(Pack pack, size_t index) noexcept;
        void store(Pack pack, size_t index, size_t count) noexcept;
    };

    template<Vector V, Packet Pack>
    void LVectorPacker<V, Pack>::store(Pack pack, size_t index) noexcept {
        using U = Traits<Pack>::ScalarType;
        if constexpr (U::isForwardDiff) {
            for (size_t i = 0; i < Pack::size(); ++i, ++index)
                (*this)[index] = pack[i];
        }
        else {
            Array<U, Pack::size()> buffer{};
            pack.store(buffer.data());
            for (size_t i = 0; i < Pack::size(); ++i, ++index)
                (*this)[index] = buffer[i];
        }
    }

    template<Vector V, Packet Pack>
    void LVectorPacker<V, Pack>::store(Pack pack, size_t index, size_t count) noexcept {
        using U = Traits<Pack>::ScalarType;
        if constexpr (U::isForwardDiff) {
            for (size_t i = 0; i < count; ++i, ++index)
                (*this)[index] = pack[i];
        }
        else {
            Array<U, Pack::size()> buffer{};
            pack.store(buffer.data());
            for (size_t i = 0; i < count; ++i, ++index)
                (*this)[index] = buffer[i];
        }
    }
}
