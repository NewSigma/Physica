/*
 * Copyright 2021-2025 Weibo He.
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
    template<class Derived> class LValueVector;
    /**
     * Reference a part of the given vector
     */
    template<Vector T, size_t Length>
    class LVectorBlock : public LValueVector<LVectorBlock<T, Length>> {
        using This = LVectorBlock<T, Length>;
        using Base = LValueVector<This>;
    public:
        using ScalarType = Base::ScalarType;
    protected:
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
    private:
        T& vec;
        size_t from;
        size_t to;
    public:
        LVectorBlock(T& vec_, size_t from_, size_t to_);
        LVectorBlock(T& vec_, size_t from_);
        LVectorBlock(const This& block) = default;
        LVectorBlock(This&&) noexcept = default;
        ~LVectorBlock() = default;
        /* Operators */
        using Base::operator=;
        LVectorBlock& operator=(const This& v) { return Base::operator=(v); }
        LVectorBlock& operator=(This&& v) noexcept { return Base::operator=(v); }
        /* Operations */
        void resize([[maybe_unused]] size_t length) const { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept;
        [[nodiscard]] PtrTy data_ptr(size_t index);
        [[nodiscard]] ConstPtrTy data_ptr(size_t index) const;
    };

    template<Vector T, size_t Length>
    LVectorBlock<T, Length>::LVectorBlock(T& vec_, size_t from_, size_t to_)
            : vec(vec_), from(from_), to(to_) {
        assert(from_ < to);
        assert(to <= vec.getLength());
        assert(Length == Dynamic || Length == getLength());
    }

    template<Vector T, size_t Length>
    LVectorBlock<T, Length>::LVectorBlock(T& vec_, size_t from_) : LVectorBlock(vec_, from_, vec_.getLength()) {}

    template<Vector T, size_t Length>
    size_t LVectorBlock<T, Length>::getLength() const noexcept {
        if constexpr (Length == Dynamic)
            return to - from;
        else
            return Length;
    }

    template<Vector T, size_t Length>
    inline auto LVectorBlock<T, Length>::data_ptr(size_t index) -> PtrTy {
        assert((index + from) < to);
        return vec.data_ptr(index + from);
    }

    template<Vector T, size_t Length>
    inline auto LVectorBlock<T, Length>::data_ptr(size_t index) const -> ConstPtrTy {
        assert((index + from) < to);
        return vec.data_ptr(index + from);
    }
}

namespace Physica {
    template<Vector T, size_t Length>
    class Traits<LVectorBlock<T, Length>> {
    public:
        using ScalarType = T::ScalarType;
        constexpr static size_t SizeAtCompile = Length;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };
}
