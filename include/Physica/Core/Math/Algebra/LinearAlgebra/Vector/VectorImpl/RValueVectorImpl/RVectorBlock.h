/*
 * Copyright 2022-2025 Weibo He.
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
    template<Vector T, size_t Length>
    class RVectorBlock : public RValueVector<RVectorBlock<T, Length>> {
        using This = RVectorBlock<T, Length>;
        using Base = RValueVector<This>;
    public:
        using ScalarType = Base::ScalarType;
    private:
        const T& vec;
        size_t from;
        size_t to;
    public:
        RVectorBlock(const T& vec_, size_t from_, size_t to_);
        RVectorBlock(const T& vec_, size_t from_);
        RVectorBlock(const This& block) = delete;
        RVectorBlock(This&&) noexcept = delete;
        ~RVectorBlock() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t index) const { return vec.calc(index + from); }
        [[nodiscard]] auto calc_value(size_t index) const { return vec.calc_value(index + from); }
        /* Getters */
        [[nodiscard]] __host__ __device__ inline size_t getLength() const noexcept;
    };

    template<Vector T, size_t Length>
    RVectorBlock<T, Length>::RVectorBlock(const T& vec_, size_t from_, size_t to_)
            : vec(vec_), from(from_), to(to_) {
        assert(from_ < to);
        assert(to <= vec.getLength());
        assert(Length == Dynamic || Length == getLength());
    }

    template<Vector T, size_t Length>
    RVectorBlock<T, Length>::RVectorBlock(const T& vec_, size_t from_) : RVectorBlock(vec_, from_, vec_.getLength()) {}

    template<Vector T, size_t Length>
    __host__ __device__ inline size_t RVectorBlock<T, Length>::getLength() const noexcept {
        if constexpr (Length == Dynamic)
            return to - from;
        else
            return Length;
    }
}

namespace Physica {
    template<Vector T, size_t Length>
    class Traits<RVectorBlock<T, Length>> {
    public:
        using ScalarType = T::ScalarType;
        constexpr static size_t SizeAtCompile = Length;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };
}
