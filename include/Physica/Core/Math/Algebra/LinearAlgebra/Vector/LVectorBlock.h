/*
 * Copyright 2021-2022 WeiBo He.
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
    class LVectorBlock;

    template<class Derived>
    class LValueVector;

    namespace Internal {
        template<class T>
        class Traits;

        template<class VectorType>
        class Traits<LVectorBlock<VectorType>> {
        public:
            using ScalarType = typename VectorType::ScalarType;
            constexpr static size_t SizeAtCompile = Dynamic;
            constexpr static size_t MaxSizeAtCompile = Dynamic;
            using PacketType = typename Internal::BestPacket<ScalarType, SizeAtCompile>::Type;
        };
    }

    template<class VectorType>
    class LVectorBlock : public LValueVector<LVectorBlock<VectorType>> {
    public:
        using Base = LValueVector<LVectorBlock<VectorType>>;
        using ScalarType = typename VectorType::ScalarType;
    private:
        VectorType& vec;
        size_t from;
        size_t to;
    public:
        LVectorBlock(LValueVector<VectorType>& vec_, size_t from_, size_t to_);
        LVectorBlock(LValueVector<VectorType>& vec_, size_t from_);
        LVectorBlock(const LVectorBlock& block) = delete;
        LVectorBlock(LVectorBlock&&) noexcept = delete;
        ~LVectorBlock() = default;
        /* Operators */
        using Base::operator=;
        LVectorBlock& operator=(const LVectorBlock& v) { Base::operator=(static_cast<const typename Base::Base&>(v)); return *this; }
        LVectorBlock& operator=(LVectorBlock&& v) noexcept { Base::operator=(static_cast<const typename Base::Base&>(v)); return *this; }
        ScalarType& operator[](size_t index) { assert((index + from) < to); return vec[index + from]; }
        const ScalarType& operator[](size_t index) const { assert((index + from) < to); return vec[index + from]; }
        /* Operations */
        void resize([[maybe_unused]] size_t length) const { assert(length == Base::getLength()); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return to - from; }
    };

    template<class VectorType>
    LVectorBlock<VectorType>::LVectorBlock(LValueVector<VectorType>& vec_, size_t from_, size_t to_)
            : vec(vec_.getDerived()), from(from_), to(to_) {
        assert(from_ < to);
        assert(to <= vec.getLength());
    }

    template<class VectorType>
    LVectorBlock<VectorType>::LVectorBlock(LValueVector<VectorType>& vec_, size_t from_) : LVectorBlock(vec_, from_, vec_.getLength()) {}
}
