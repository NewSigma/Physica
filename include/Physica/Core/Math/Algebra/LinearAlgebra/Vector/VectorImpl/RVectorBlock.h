/*
 * Copyright 2022 WeiBo He.
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
    using Utils::Dynamic;
    /**
     * Reference a part of the given vector
     */
    template<class T>
    class RVectorBlock;

    template<class Derived>
    class RValueVector;

    namespace Internal {
        template<class T>
        class Traits;

        template<class VectorType>
        class Traits<RVectorBlock<VectorType>> {
        public:
            using ScalarType = typename VectorType::ScalarType;
            constexpr static size_t SizeAtCompile = Dynamic;
            constexpr static size_t MaxSizeAtCompile = Dynamic;
            using PacketType = typename Internal::BestPacket<ScalarType, SizeAtCompile>::Type;
        };
    }

    template<class VectorType>
    class RVectorBlock : public RValueVector<RVectorBlock<VectorType>> {
    public:
        using Base = RValueVector<RVectorBlock<VectorType>>;
        using ScalarType = typename VectorType::ScalarType;
    private:
        const VectorType& vec;
        size_t from;
        size_t to;
    public:
        RVectorBlock(const RValueVector<VectorType>& vec_, size_t from_, size_t to_);
        RVectorBlock(const RValueVector<VectorType>& vec_, size_t from_);
        RVectorBlock(const RVectorBlock& block) = delete;
        RVectorBlock(RVectorBlock&&) noexcept = delete;
        ~RVectorBlock() = default;
        /* Operators */
        [[nodiscard]] ScalarType calc(size_t index) const { return vec.calc(index + from); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return to - from; }
    };

    template<class VectorType>
    RVectorBlock<VectorType>::RVectorBlock(const RValueVector<VectorType>& vec_, size_t from_, size_t to_)
            : vec(vec_.getDerived()), from(from_), to(to_) {
        assert(from_ < to);
        assert(to <= vec.getLength());
    }

    template<class VectorType>
    RVectorBlock<VectorType>::RVectorBlock(const RValueVector<VectorType>& vec_, size_t from_) : RVectorBlock(vec_, from_, vec_.getLength()) {}
}
