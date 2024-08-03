/*
 * Copyright 2021-2024 Weibo He.
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
    template<class Derived> class LValueVector;
    /**
     * Reference a part of the given vector
     */
    template<class VectorType>
    class LVectorBlock : public LValueVector<LVectorBlock<VectorType>> {
        using This = LVectorBlock<VectorType>;
    public:
        using Base = LValueVector<This>;
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
        LVectorBlock& operator=(const LVectorBlock& v) { Base::operator=(static_cast<const RValueVector<This>&>(v)); return *this; }
        LVectorBlock& operator=(LVectorBlock&& v) noexcept { Base::operator=(static_cast<const RValueVector<This>&>(v)); return *this; }
        /* Operations */
        void resize([[maybe_unused]] size_t length) const { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return to - from; }
        [[nodiscard]] __host__ __device__ ScalarType* data_ptr(size_t index);
        [[nodiscard]] __host__ __device__ const ScalarType* data_ptr(size_t index) const;
    };

    template<class VectorType>
    LVectorBlock<VectorType>::LVectorBlock(LValueVector<VectorType>& vec_, size_t from_, size_t to_)
            : vec(vec_.getDerived()), from(from_), to(to_) {
        assert(from_ < to);
        assert(to <= vec.getLength());
    }

    template<class VectorType>
    LVectorBlock<VectorType>::LVectorBlock(LValueVector<VectorType>& vec_, size_t from_) : LVectorBlock(vec_, from_, vec_.getLength()) {}

    template<class VectorType>
    __host__ __device__ inline typename LVectorBlock<VectorType>::ScalarType*
    LVectorBlock<VectorType>::data_ptr(size_t index) {
        assert((index + from) < to);
        return vec.data_ptr(index + from);
    }

    template<class VectorType>
    __host__ __device__ inline const typename LVectorBlock<VectorType>::ScalarType*
    LVectorBlock<VectorType>::data_ptr(size_t index) const {
        assert((index + from) < to);
        return vec.data_ptr(index + from);
    }
}

namespace Physica {
    using namespace Core;

    template<class VectorType>
    class Traits<LVectorBlock<VectorType>> {
    public:
        using ScalarType = typename VectorType::ScalarType;
        constexpr static size_t SizeAtCompile = Dynamic;
        constexpr static size_t MaxSizeAtCompile = Dynamic;

        constexpr static bool FastAssign = false;
    };
}
