/*
 * Copyright 2026 Weibo He.
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

#include "../RValueVector.h"

namespace Physica {
    /**
     * \class RVectorPacker is base class of all packers. Packer is a range for SIMD packets.
     *
     * It is basically a loop invariant code motion. We do this because compilers struggle to prove noalias.
     */
    template<Vector V, Packet Pack>
    class RVectorPacker : public std::ranges::view_base {
        using This = RVectorPacker<V, Pack>;
        using Base = CRTPBase<This>;
    private:
        class Iterator : public std::forward_iterator_tag {
            using This = Iterator;

            const RVectorPacker<V, Pack>* packer;
            size_t index{};
        public:
            using difference_type = std::ptrdiff_t;
            using value_type = const Pack;
            using reference = value_type;
            using const_reference = value_type;
        public:
            Iterator() = default;
            Iterator(const RVectorPacker<V, Pack>* packer, size_t index);
            Iterator(const This&) = default;
            Iterator(This&&) noexcept = default;
            ~Iterator() = default;
            /* Operators */
            This& operator=(const This&) = default;
            This& operator=(This&&) noexcept = default;
            [[nodiscard]] bool operator==(const This& other) const noexcept;
            [[nodiscard]] auto operator<=>(const This& other) const noexcept;
            This& operator++() noexcept;
            This operator++(int) noexcept;
            [[nodiscard]] reference operator*() const noexcept;
        };

        const V* vec;
    public:
        RVectorPacker(const V& vec) : vec(&vec) {}
        RVectorPacker(const This&) = default;
        RVectorPacker(This&&) noexcept = default;
        ~RVectorPacker() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        [[nodiscard]] Pack load(size_t index) const noexcept;
        [[nodiscard]] Pack load(size_t index, size_t count) const noexcept;
        [[nodiscard]] Pack tail() const noexcept;

        [[nodiscard]] auto begin() noexcept { return cbegin(); }
        [[nodiscard]] auto begin() const noexcept { return cbegin(); }
        [[nodiscard]] auto cbegin() const noexcept { return Iterator(this, 0); }
        [[nodiscard]] auto end() noexcept { return cend(); }
        [[nodiscard]] auto end() const noexcept { return cend(); }
        [[nodiscard]] auto cend() const noexcept { return Iterator(this, size() * Pack::size()); }
        /* Getters */
        [[nodiscard]] const auto& vector() const noexcept { return *vec; }
        [[nodiscard]] constexpr size_t size() const noexcept;
        [[nodiscard]] size_t tailing() const noexcept;
        /* Static members */
        [[nodiscard]] static Pack zeros() noexcept { return Pack::zeros(); }
    };

    template<Vector V, Packet Pack>
    RVectorPacker<V, Pack>::Iterator::Iterator(const RVectorPacker<V, Pack>* packer, size_t index) : packer(packer), index(index) {}

    template<Vector V, Packet Pack>
    bool RVectorPacker<V, Pack>::Iterator::operator==(const This& other) const noexcept {
        return index == other.index;
    }

    template<Vector V, Packet Pack>
    auto RVectorPacker<V, Pack>::Iterator::operator<=>(const This& other) const noexcept {
        return index <=> other.index;
    }

    template<Vector V, Packet Pack>
    auto RVectorPacker<V, Pack>::Iterator::operator++() noexcept -> This& {
        index += Pack::size();
        return *this;
    }

    template<Vector V, Packet Pack>
    auto RVectorPacker<V, Pack>::Iterator::operator++(int) noexcept -> This {
        return std::exchange(*this, This(packer, index + Pack::size()));
    }

    template<Vector V, Packet Pack>
    auto RVectorPacker<V, Pack>::Iterator::operator*() const noexcept -> reference {
        return packer->load(index);
    }

    template<Vector V, Packet Pack>
    auto RVectorPacker<V, Pack>::load(size_t index) const noexcept -> Pack {
        using U = Traits<Pack>::ScalarType;
        assert(index + Pack::size() <= vec->getLength() && "[Error]: Index out of range");
        if constexpr (Diffable<U>) {
            using ValuePacket = Pack::ValueType;
            if constexpr (V::isForwardDiff) {
                using GradPacket = Pack::GradType;
                auto values = vec->values().template packet<ValuePacket>(index);
                auto grads = vec->grads().template packet<GradPacket>(index);
                return Pack(std::move(values), std::move(grads));
            }
            else {
                using Uv = U::ValueType;
                Array<Uv, Pack::size()> values{};
                for (size_t i = 0; i < Pack::size(); ++i, ++index)
                    values[i] = Uv(vec->calc(index));
                ValuePacket packet{};
                packet.load(values.data());
                return Pack(std::move(packet));
            }
        }
        else {
            Array<U, Pack::size()> buffer{};
            for (size_t i = 0; i < Pack::size(); ++i, ++index)
                buffer[i] = U(vec->calc(index));
            Pack packet{};
            packet.load(buffer.data());
            return packet;
        }
    }

    template<Vector V, Packet Pack>
    auto RVectorPacker<V, Pack>::load(size_t index, size_t count) const noexcept -> Pack {
        using U = Traits<Pack>::ScalarType;
        assert(index + count <= vec->getLength() && "[Error]: Index out of range");
        if constexpr (Diffable<U>) {
            using ValuePacket = Pack::ValueType;
            if constexpr (V::isForwardDiff) {
                using GradPacket = Pack::GradType;
                auto values = vec->values().template packet<ValuePacket>(index, count);
                auto grads = vec->grads().template packet<GradPacket>(index, count);
                return Pack(std::move(values), std::move(grads));
            }
            else {
                using Uv = U::ValueType;
                Array<Uv, Pack::size()> values{};
                for (size_t i = 0; i < Pack::size(); ++i, ++index)
                    values[i] = i < count ? Uv(vec->calc(index)) : Uv(0);
                ValuePacket packet{};
                packet.load(values.data());
                return Pack(std::move(packet));
            }
        }
        else {
            Array<U, Pack::size()> buffer{};
            for (size_t i = 0; i < Pack::size(); ++i, ++index)
                buffer[i] = i < count ? U(vec->calc(index)) : U(0);
            Pack packet{};
            packet.load(buffer.data());
            return packet;
        }
    }

    template<Vector V, Packet Pack>
    auto RVectorPacker<V, Pack>::tail() const noexcept -> Pack {
        size_t index = size() * Pack::size();
        return load(index, vec->getLength() - index);
    }

    template<Vector V, Packet Pack>
    constexpr size_t RVectorPacker<V, Pack>::size() const noexcept {
        if constexpr (V::SizeAtCompile != Dynamic)
            return V::SizeAtCompile / Pack::size();
        else
            return vec->getLength() / Pack::size();
    }

    template<Vector V, Packet Pack>
    size_t RVectorPacker<V, Pack>::tailing() const noexcept {
        return vec->getLength() % Pack::size();
    }
}
