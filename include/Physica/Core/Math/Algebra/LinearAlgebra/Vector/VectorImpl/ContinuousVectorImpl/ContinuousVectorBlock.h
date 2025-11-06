/*
 * Copyright 2023-2025 Weibo He.
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
    template<class Derived> class ContinuousVector;

    template<Vector T, size_t Length>
    class ContinuousVectorBlock : public ContinuousVector<ContinuousVectorBlock<T, Length>> {
        using This = ContinuousVectorBlock<T, Length>;
        using Base = ContinuousVector<This>;
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
        ContinuousVectorBlock(const This& block) = default;
        ContinuousVectorBlock(This&&) noexcept = default;
        ~ContinuousVectorBlock() = default;
        /* Operators */
        using Base::operator=;
        This& operator=(const This& v);
        This& operator=(This&& v) noexcept;
        [[nodiscard]] RefTy operator[](size_t index);
        [[nodiscard]] ConstRefTy operator[](size_t index) const;
        /* Operations */
        template<Packet Pack> [[nodiscard]] Pack packet(size_t index) const;
        template<Packet Pack> [[nodiscard]] Pack packetPartial(size_t index, size_t count) const;
        void writePacket(size_t index, Packet auto packet);
        void writePacketPartial(size_t index, size_t count, Packet auto packet);

        using Base::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept;
        [[nodiscard]] PtrTy data_ptr(size_t index) noexcept { return vec.data() + from + index; }
        [[nodiscard]] ConstPtrTy data_ptr(size_t index) const noexcept { return vec.data() + from + index; }
    };

    template<Vector T, size_t Length>
    ContinuousVectorBlock<T, Length>::ContinuousVectorBlock(ContinuousVector<T>& vec_, size_t from_, size_t to_)
            : vec(vec_.getDerived()), from(from_), to(to_) {
        assert(from_ < to);
        assert(to <= vec.getLength());
        assert(Length == Dynamic || Length == getLength());
    }

    template<Vector T, size_t Length>
    ContinuousVectorBlock<T, Length>::ContinuousVectorBlock(ContinuousVector<T>& vec_, size_t from_)
            : ContinuousVectorBlock(vec_, from_, vec_.getLength()) {}

    template<Vector T, size_t Length>
    auto ContinuousVectorBlock<T, Length>::operator=(const This& v) -> This& {
        Base::operator=(v);
        return *this;
    }

    template<Vector T, size_t Length>
    auto ContinuousVectorBlock<T, Length>::operator=(This&& v) noexcept -> This& {
        Base::operator=(v);
        return *this;
    }

    template<Vector T, size_t Length>
    auto ContinuousVectorBlock<T, Length>::operator[](size_t index) -> RefTy {
        assert((index + from) < to);
        return vec[index + from];
    }

    template<Vector T, size_t Length>
    auto ContinuousVectorBlock<T, Length>::operator[](size_t index) const -> ConstRefTy {
        assert((index + from) < to);
        return vec[index + from];
    }

    template<Vector T, size_t Length>
    template<Packet Pack>
    Pack ContinuousVectorBlock<T, Length>::packet(size_t index) const {
        return vec.template packet<Pack>(from + index);
    }

    template<Vector T, size_t Length>
    template<Packet Pack>
    Pack ContinuousVectorBlock<T, Length>::packetPartial(size_t index, size_t count) const {
        return vec.template packetPartial<Pack>(from + index, count);
    }

    template<Vector T, size_t Length>
    void ContinuousVectorBlock<T, Length>::writePacket(size_t index, const Packet auto packet) {
        return vec.writePacket(from + index, packet);
    }

    template<Vector T, size_t Length>
    void ContinuousVectorBlock<T, Length>::writePacketPartial(size_t index, size_t count, const Packet auto packet) {
        return vec.writePacketPartial(from + index, count, packet);
    }

    template<Vector T, size_t Length>
    size_t ContinuousVectorBlock<T, Length>::getLength() const noexcept {
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

        using ElemType = ScalarType;
    };
}
