/*
 * Copyright 2023-2024 Weibo He.
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
    template<class Derived> class ContinuousVector;

    template<Vector T, size_t Length>
    class ContinuousVectorBlock : public ContinuousVector<ContinuousVectorBlock<T, Length>> {
        using This = ContinuousVectorBlock<T, Length>;
    public:
        using Base = ContinuousVector<This>;
        using ScalarType = T::ScalarType;
    protected:
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
        using typename Base::RefTy;
        using typename Base::ConstRefTy;
    private:
        T& vec;
        size_t from;
        size_t to;
    public:
        ContinuousVectorBlock(ContinuousVector<T>& vec_, size_t from_, size_t to_);
        ContinuousVectorBlock(ContinuousVector<T>& vec_, size_t from_);
        ContinuousVectorBlock(const ContinuousVectorBlock& block) = delete;
        ContinuousVectorBlock(ContinuousVectorBlock&&) noexcept = delete;
        ~ContinuousVectorBlock() = default;
        /* Operators */
        using Base::operator=;
        inline ContinuousVectorBlock& operator=(const ContinuousVectorBlock& v);
        inline ContinuousVectorBlock& operator=(ContinuousVectorBlock&& v) noexcept;
        [[nodiscard]] RefTy operator[](size_t index);
        [[nodiscard]] ConstRefTy operator[](size_t index) const;
        /* Operations */
        template<class AnyPacket> [[nodiscard]] inline AnyPacket packet(size_t index) const;
        template<class AnyPacket> [[nodiscard]] inline AnyPacket packetPartial(size_t index, size_t count) const;
        template<class AnyPacket> inline void writePacket(size_t index, const AnyPacket packet);
        template<class AnyPacket> inline void writePacketPartial(size_t index, size_t count, const AnyPacket packet);
        void resize([[maybe_unused]] size_t length) const { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ inline size_t getLength() const noexcept;
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(size_t index) { return vec.data() + from + index; }
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(size_t index) const { return vec.data() + from + index; }
    };

    template<Vector T, size_t Length>
    ContinuousVectorBlock<T, Length>::ContinuousVectorBlock(ContinuousVector<T>& vec_, size_t from_, size_t to_)
            : vec(vec_.getDerived()), from(from_), to(to_) {
        assert(from_ < to);
        assert(to <= vec.getLength());
        assert(Length == Dynamic || Length == getLength());
    }

    template<Vector T, size_t Length>
    ContinuousVectorBlock<T, Length>::ContinuousVectorBlock(
        ContinuousVector<T>& vec_, size_t from_) : ContinuousVectorBlock(vec_, from_, vec_.getLength()) {}

    template<Vector T, size_t Length>
    inline ContinuousVectorBlock<T, Length>&
    ContinuousVectorBlock<T, Length>::operator=(const ContinuousVectorBlock<T, Length>& v) {
        v.assignTo(*this);
        return *this;
    }
    
    template<Vector T, size_t Length>
    inline ContinuousVectorBlock<T, Length>&
    ContinuousVectorBlock<T, Length>::operator=(ContinuousVectorBlock<T, Length>&& v) noexcept {
        v.assignTo(*this);
        return *this;
    }

    template<Vector T, size_t Length>
    inline ContinuousVectorBlock<T, Length>::RefTy
    ContinuousVectorBlock<T, Length>::operator[](size_t index) {
        assert((index + from) < to);
        return vec[index + from];
    }

    template<Vector T, size_t Length>
    inline ContinuousVectorBlock<T, Length>::ConstRefTy
    ContinuousVectorBlock<T, Length>::operator[](size_t index) const {
        assert((index + from) < to);
        return vec[index + from];
    }

    template<Vector T, size_t Length>
    template<class AnyPacket>
    inline AnyPacket ContinuousVectorBlock<T, Length>::packet(size_t index) const {
        return vec.template packet<AnyPacket>(from + index);
    }

    template<Vector T, size_t Length>
    template<class AnyPacket>
    inline AnyPacket ContinuousVectorBlock<T, Length>::packetPartial(size_t index, size_t count) const {
        return vec.template packetPartial<AnyPacket>(from + index, count);
    }

    template<Vector T, size_t Length>
    template<class AnyPacket>
    inline void ContinuousVectorBlock<T, Length>::writePacket(size_t index, const AnyPacket packet) {
        return vec.template writePacket<AnyPacket>(from + index, packet);
    }

    template<Vector T, size_t Length>
    template<class AnyPacket>
    inline void ContinuousVectorBlock<T, Length>::writePacketPartial(size_t index, size_t count, const AnyPacket packet) {
        return vec.template writePacketPartial<AnyPacket>(from + index, count, packet);
    }

    template<Vector T, size_t Length>
    __host__ __device__ inline size_t ContinuousVectorBlock<T, Length>::getLength() const noexcept {
        if constexpr (Length == Dynamic)
            return to - from;
        else
            return Length;
    }
}

namespace Physica {
    template<Vector T, size_t Length>
    class Traits<ContinuousVectorBlock<T, Length>> {
    public:
        using ScalarType = T::ScalarType;
        constexpr static size_t SizeAtCompile = Length;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = true;
    };
}
