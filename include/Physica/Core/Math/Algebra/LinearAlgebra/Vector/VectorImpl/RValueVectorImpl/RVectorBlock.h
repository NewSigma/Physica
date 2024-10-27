/*
 * Copyright 2022-2024 Weibo He.
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

namespace Physica::Core {
    template<class Derived> class RValueVector;
    /**
     * Reference a part of the given vector
     */
    template<class VectorType, size_t Length>
    class RVectorBlock : public RValueVector<RVectorBlock<VectorType, Length>> {
        using This = RVectorBlock<VectorType, Length>;
        using Base = RValueVector<This>;
    public:
        using ScalarType = typename Base::ScalarType;
    private:
        const VectorType& vec;
        size_t from;
        size_t to;
    public:
        RVectorBlock(const RValueVector<VectorType>& vec_, size_t from_, size_t to_);
        RVectorBlock(const RValueVector<VectorType>& vec_, size_t from_);
        RVectorBlock(const This& block) = delete;
        RVectorBlock(This&&) noexcept = delete;
        ~RVectorBlock() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t index) const { return vec.calc(index + from); }
        /* Getters */
        [[nodiscard]] __host__ __device__ inline size_t getLength() const noexcept;
    };

    template<class VectorType, size_t Length>
    RVectorBlock<VectorType, Length>::RVectorBlock(const RValueVector<VectorType>& vec_, size_t from_, size_t to_)
            : vec(vec_.getDerived()), from(from_), to(to_) {
        assert(from_ < to);
        assert(to <= vec.getLength());
        assert(Length == Dynamic || Length == getLength());
    }

    template<class VectorType, size_t Length>
    RVectorBlock<VectorType, Length>::RVectorBlock(const RValueVector<VectorType>& vec_, size_t from_) : RVectorBlock(vec_, from_, vec_.getLength()) {}

    template<class VectorType, size_t Length>
    __host__ __device__ inline size_t RVectorBlock<VectorType, Length>::getLength() const noexcept {
        if constexpr (Length == Dynamic)
            return to - from;
        else
            return Length;
    }
}

namespace Physica {
    template<class VectorType, size_t Length>
    class Traits<Core::RVectorBlock<VectorType, Length>> {
    public:
        using ScalarType = typename VectorType::ScalarType;
        constexpr static size_t SizeAtCompile = Length;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };
}
