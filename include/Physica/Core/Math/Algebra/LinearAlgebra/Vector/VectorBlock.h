/*
 * Copyright 2021 WeiBo He.
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
    /**
     * Reference a part of the given vector
     */
    template<class T>
    class VectorBlock;

    template<class Derived>
    class LValueVector;

    namespace Internal {
        template<class T>
        class Traits;

        template<class VectorType>
        class Traits<VectorBlock<VectorType>> {
        public:
            using ScalarType = typename VectorType::ScalarType;
        };
    }

    template<class VectorType>
    class VectorBlock : public LValueVector<VectorBlock<VectorType>> {
    public:
        using ScalarType = typename VectorType::ScalarType;
    private:
        VectorType& vec;
        size_t from;
        size_t to;
    public:
        VectorBlock(LValueVector<VectorType>& vec_, size_t from_, size_t to_);
        VectorBlock(LValueVector<VectorType>& vec_, size_t from_);
        VectorBlock(const VectorBlock& block);
        VectorBlock(VectorBlock&&) noexcept = delete;
        ~VectorBlock() = default;
        /* Operators */
        VectorBlock& operator=(const VectorBlock&) = delete;
        VectorBlock& operator=(VectorBlock&&) noexcept = delete;
        template<class OtherVector>
        VectorBlock& operator=(const RValueVector<OtherVector>& v);
        ScalarType& operator[](size_t index) { assert((index + from) < to); return vec[index + from]; }
        const ScalarType& operator[](size_t index) const { assert((index + from) < to); return vec[index + from]; }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return to - from; }
    };

    template<class VectorType>
    VectorBlock<VectorType>::VectorBlock(LValueVector<VectorType>& vec_, size_t from_, size_t to_)
            : vec(vec_.getDerived()), from(from_), to(to_) {
        assert(from_ < to_);
        assert(to_ <= vec.getLength());
    }

    template<class VectorType>
    VectorBlock<VectorType>::VectorBlock(LValueVector<VectorType>& vec_, size_t from_) : VectorBlock(vec_, from_, vec_.getLength()) {}

    template<class VectorType>
    template<class OtherVector>
    VectorBlock<VectorType>& VectorBlock<VectorType>::operator=(const RValueVector<OtherVector>& v) {
        assert(getLength() == v.getLength());
        const size_t length = getLength();
        for (size_t i = 0; i < length; ++i)
            this->operator[](i) = v.calc(i);
        return *this;
    }

    template<class VectorType>
    VectorBlock<VectorType>::VectorBlock(const VectorBlock& block) : vec(block.vec), from(block.from), to(block.to) {}
}
