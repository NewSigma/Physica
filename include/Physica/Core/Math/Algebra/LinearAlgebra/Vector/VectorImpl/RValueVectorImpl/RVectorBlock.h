/*
 * Copyright 2022-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

namespace Physica {
    /**
     * Reference a part of the given vector
     */
    template<Vector V, size_t Length>
    class RVectorBlock : public RValueVector<RVectorBlock<V, Length>> {
        using This = RVectorBlock<V, Length>;
        using Base = RValueVector<This>;
    protected:
        using typename Base::T;
    private:
        V& vec;
        size_t from;
        size_t to;
    public:
        RVectorBlock(V& vec_, size_t from_, size_t to_);
        RVectorBlock(V& vec_, size_t from_);
        RVectorBlock(const This& block) = default;
        RVectorBlock(This&&) noexcept = default;
        ~RVectorBlock() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t index) const { return vec.calc(index + from); }
        [[nodiscard]] auto calc_value(size_t index) const { return vec.calc_value(index + from); }

        template<size_t Length_ = Dynamic>
        [[nodiscard]] auto head(this auto&&, size_t to = Length_) noexcept;
        template<size_t Length_ = Dynamic>
        [[nodiscard]] auto tail(this auto&&, size_t from) noexcept;
        template<size_t Length_ = Dynamic>
        [[nodiscard]] auto segment(this auto&&, size_t from, size_t to) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept;
    };

    template<Vector V, size_t Length>
    RVectorBlock<V, Length>::RVectorBlock(V& vec_, size_t from_, size_t to_) : vec(vec_), from(from_), to(to_) {
        assert(from_ < to);
        assert(to <= vec.getLength());
        assert(Length == Dynamic || Length == getLength());
    }

    template<Vector V, size_t Length>
    RVectorBlock<V, Length>::RVectorBlock(V& vec_, size_t from_) : RVectorBlock(vec_, from_, vec_.getLength()) {}

    template<Vector V, size_t Length>
    template<size_t Length_>
    auto RVectorBlock<V, Length>::head(this auto&& self, size_t to) noexcept {
        return RVectorBlock<V, Length_>(self.vec, self.from, self.from + to);
    }

    template<Vector V, size_t Length>
    template<size_t Length_>
    auto RVectorBlock<V, Length>::tail(this auto&& self, size_t from) noexcept {
        return RVectorBlock<V, Length_>(self.vec, self.from + from, self.to);
    }

    template<Vector V, size_t Length>
    template<size_t Length_>
    auto RVectorBlock<V, Length>::segment(this auto&& self, size_t from, size_t to) noexcept {
        return RVectorBlock<V, Length_>(self.vec, self.from + from, self.from + to);
    }

    template<Vector V, size_t Length>
    size_t RVectorBlock<V, Length>::getLength() const noexcept {
        if constexpr (Length == Dynamic)
            return to - from;
        else
            return Length;
    }
}

namespace Physica {
    template<Vector V, size_t Length>
    class Traits<RVectorBlock<V, Length>> {
        static_assert(std::is_object<V>::value, "[Errpr]: Invalid type");
    public:
        using ScalarType = V::ScalarType;
        constexpr static size_t SizeAtCompile = Length;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };
}
