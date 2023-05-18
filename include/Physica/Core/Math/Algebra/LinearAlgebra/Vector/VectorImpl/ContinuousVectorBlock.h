/*
 * Copyright 2023 WeiBo He.
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
    template<class VectorType, size_t Length> class ContinuousVectorBlock;

    namespace Internal {
        template<class T> class Traits;

        template<class VectorType, size_t Length>
        class Traits<ContinuousVectorBlock<VectorType, Length>> {
        public:
            using ScalarType = typename VectorType::ScalarType;
            constexpr static size_t SizeAtCompile = Length;
            constexpr static size_t MaxSizeAtCompile = Length;
            using PacketType = typename Internal::BestPacket<ScalarType, SizeAtCompile>::Type;
        };
    }

    template<class VectorType, size_t Length>
    class ContinuousVectorBlock : public ContinuousVector<ContinuousVectorBlock<VectorType, Length>> {
        using This = ContinuousVectorBlock<VectorType, Length>;
    public:
        using Base = ContinuousVector<This>;
        using ScalarType = typename VectorType::ScalarType;
    private:
        VectorType& vec;
        size_t from;
        size_t to;
    public:
        ContinuousVectorBlock(ContinuousVector<VectorType>& vec_, size_t from_, size_t to_);
        ContinuousVectorBlock(ContinuousVector<VectorType>& vec_, size_t from_);
        ContinuousVectorBlock(const ContinuousVectorBlock& block) = delete;
        ContinuousVectorBlock(ContinuousVectorBlock&&) noexcept = delete;
        ~ContinuousVectorBlock() = default;
        /* Operators */
        using Base::operator=;
        inline ContinuousVectorBlock& operator=(const ContinuousVectorBlock& v);
        inline ContinuousVectorBlock& operator=(ContinuousVectorBlock&& v) noexcept;
        [[nodiscard]] ScalarType& operator[](size_t index);
        [[nodiscard]] const ScalarType& operator[](size_t index) const;
        /* Operations */
        template<class PacketType> [[nodiscard]] inline PacketType packet(size_t index) const;
        template<class PacketType> [[nodiscard]] inline PacketType packetPartial(size_t index, size_t count) const;
        template<class PacketType> inline void writePacket(size_t index, const PacketType packet);
        template<class PacketType> inline void writePacketPartial(size_t index, size_t count, const PacketType packet);
        void resize([[maybe_unused]] size_t length) const { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] inline size_t getLength() const noexcept;
    };

    template<class VectorType, size_t Length>
    ContinuousVectorBlock<VectorType, Length>::ContinuousVectorBlock(ContinuousVector<VectorType>& vec_, size_t from_, size_t to_)
            : vec(vec_.getDerived()), from(from_), to(to_) {
        assert(from_ < to);
        assert(to <= vec.getLength());
        assert(Length == Dynamic || Length == getLength());
    }

    template<class VectorType, size_t Length>
    ContinuousVectorBlock<VectorType, Length>::ContinuousVectorBlock(
        ContinuousVector<VectorType>& vec_, size_t from_) : ContinuousVectorBlock(vec_, from_, vec_.getLength()) {}

    template<class VectorType, size_t Length>
    inline ContinuousVectorBlock<VectorType, Length>&
    ContinuousVectorBlock<VectorType, Length>::operator=(const ContinuousVectorBlock<VectorType, Length>& v) {
        Base::operator=(v);
        return *this;
    }
    
    template<class VectorType, size_t Length>
    inline ContinuousVectorBlock<VectorType, Length>&
    ContinuousVectorBlock<VectorType, Length>::operator=(ContinuousVectorBlock<VectorType, Length>&& v) noexcept {
        Base::operator=(std::move(v));
        return *this;
    }

    template<class VectorType, size_t Length>
    inline typename ContinuousVectorBlock<VectorType, Length>::ScalarType&
    ContinuousVectorBlock<VectorType, Length>::operator[](size_t index) {
        assert((index + from) < to);
        return vec[index + from];
    }

    template<class VectorType, size_t Length>
    inline const typename ContinuousVectorBlock<VectorType, Length>::ScalarType&
    ContinuousVectorBlock<VectorType, Length>::operator[](size_t index) const {
        assert((index + from) < to);
        return vec[index + from];
    }

    template<class VectorType, size_t Length>
    template<class PacketType>
    inline PacketType ContinuousVectorBlock<VectorType, Length>::packet(size_t index) const {
        return vec.template packet<PacketType>(from + index);
    }

    template<class VectorType, size_t Length>
    template<class PacketType>
    inline PacketType ContinuousVectorBlock<VectorType, Length>::packetPartial(size_t index, size_t count) const {
        return vec.template packetPartial<PacketType>(from + index, count);
    }

    template<class VectorType, size_t Length>
    template<class PacketType>
    inline void ContinuousVectorBlock<VectorType, Length>::writePacket(size_t index, const PacketType packet) {
        return vec.template writePacket<PacketType>(from + index, packet);
    }

    template<class VectorType, size_t Length>
    template<class PacketType>
    inline void ContinuousVectorBlock<VectorType, Length>::writePacketPartial(size_t index, size_t count, const PacketType packet) {
        return vec.template writePacketPartial<PacketType>(from + index, count, packet);
    }

    template<class VectorType, size_t Length>
    inline size_t ContinuousVectorBlock<VectorType, Length>::getLength() const noexcept {
        if constexpr (Length == Dynamic)
            return to - from;
        return Length;
    }
}
