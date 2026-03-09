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

#include "../CompactVector.h"

namespace Physica {
    template<Vector V>
    class CompactVectorView final : public LVectorView<V> {
        using This = CompactVectorView<V>;
        using Base = LVectorView<V>;

        class Iterator;
    public:
        constexpr CompactVectorView(V& vec) noexcept;
        constexpr CompactVectorView(const This&) = default;
        constexpr CompactVectorView(This&&) noexcept = default;
        constexpr ~CompactVectorView() = default;
        /* Operators */
        constexpr This& operator=(const This&) = default;
        constexpr This& operator=(This&&) noexcept = default;
        /* Operations */
        [[nodiscard]] constexpr auto begin(this auto&&) noexcept;
        [[nodiscard]] constexpr auto end(this auto&&) noexcept;
    };

    template<Vector V>
    constexpr CompactVectorView<V>::CompactVectorView(V& vec) noexcept : Base(vec) {}

    template<Vector V>
    constexpr auto CompactVectorView<V>::begin(this auto&& self) noexcept {
        return Iterator(self.vector().data());
    }

    template<Vector V>
    constexpr auto CompactVectorView<V>::end(this auto&& self) noexcept {
        return Iterator(self.vector().data() + self.size());
    }
}

namespace Physica {
    template<Vector V>
    class CompactVectorView<V>::Iterator {
        using This = Iterator;
    public:
        using iterator_concept = std::contiguous_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = V::ScalarType;
        using reference = value_type::RefTy;
        using const_reference = value_type::ConstRefTy;
        using pointer = std::conditional<std::is_const_v<V>, typename value_type::ConstPtrTy, typename value_type::PtrTy>::type;
    private:
        pointer data;
    public:
        constexpr Iterator() = default;
        constexpr Iterator(pointer data_) noexcept;
        constexpr Iterator(const This&) = default;
        constexpr Iterator(This&&) noexcept = default;
        constexpr ~Iterator() = default;
        /* Operators */
        constexpr This& operator=(const This&) = default;
        constexpr This& operator=(This&&) noexcept = default;
        constexpr This& operator++() noexcept;
        constexpr This& operator--() noexcept;
        constexpr This& operator+=(difference_type n) noexcept;
        constexpr This& operator-=(difference_type n) noexcept;
        [[nodiscard]] constexpr This operator++(int) noexcept;
        [[nodiscard]] constexpr This operator--(int) noexcept;
        [[nodiscard]] constexpr decltype(auto) operator*() const noexcept;
        [[nodiscard]] constexpr decltype(auto) operator[](difference_type n) const noexcept;
        [[nodiscard]] constexpr bool operator==(const This& other) const noexcept;
        [[nodiscard]] constexpr auto operator<=>(const This& other) const noexcept;
        [[nodiscard]] constexpr This operator+(difference_type n) const noexcept;
        [[nodiscard]] constexpr This operator-(difference_type n) const noexcept;
        [[nodiscard]] constexpr difference_type operator-(const This& other) const noexcept;
        /* Operations */
        /* Operations */
        template<Packet Pack>
        [[nodiscard]] Pack load() const noexcept;
        template<Packet Pack>
        [[nodiscard]] Pack load(size_t count) const noexcept;
        void store(Packet auto pack) noexcept;
        void store(Packet auto pack, size_t count) noexcept;
        /* Friends */
        friend constexpr This operator+(difference_type n, const This& ite) noexcept { return ite + n; }
    };

    template<Vector V>
    constexpr CompactVectorView<V>::Iterator::Iterator(pointer data_) noexcept : data(data_) {}

    template<Vector V>
    constexpr auto CompactVectorView<V>::Iterator::operator++() noexcept -> This& {
        data += 1;
        return *this;
    }

    template<Vector V>
    constexpr auto CompactVectorView<V>::Iterator::operator--() noexcept -> This& {
        data -= 1;
        return *this;
    }

    template<Vector V>
    constexpr auto CompactVectorView<V>::Iterator::operator+=(difference_type n) noexcept -> This& {
        data += n;
        return *this;
    }

    template<Vector V>
    constexpr auto CompactVectorView<V>::Iterator::operator-=(difference_type n) noexcept -> This& {
        data -= n;
        return *this;
    }

    template<Vector V>
    constexpr auto CompactVectorView<V>::Iterator::operator++(int) noexcept -> This {
        return std::exchange(*this, This(data + 1));
    }

    template<Vector V>
    constexpr auto CompactVectorView<V>::Iterator::operator--(int) noexcept -> This {
        return std::exchange(*this, This(data - 1));
    }

    template<Vector V>
    constexpr decltype(auto) CompactVectorView<V>::Iterator::operator*() const noexcept {
        return *data;
    }

    template<Vector V>
    constexpr decltype(auto) CompactVectorView<V>::Iterator::operator[](difference_type n) const noexcept {
        return data[n];
    }

    template<Vector V>
    constexpr bool CompactVectorView<V>::Iterator::operator==(const This& other) const noexcept {
        return data == other.data;
    }

    template<Vector V>
    constexpr auto CompactVectorView<V>::Iterator::operator<=>(const This& other) const noexcept {
        return data <=> other.data;
    }

    template<Vector V>
    constexpr auto CompactVectorView<V>::Iterator::operator+(difference_type n) const noexcept -> This {
        return This(data + n);
    }

    template<Vector V>
    constexpr auto CompactVectorView<V>::Iterator::operator-(difference_type n) const noexcept -> This {
        return This(data - n);
    }

    template<Vector V>
    constexpr auto CompactVectorView<V>::Iterator::operator-(const This& other) const noexcept -> difference_type {
        return data - other.data;
    }

    template<Vector V>
    template<Packet Pack>
    auto CompactVectorView<V>::Iterator::load() const noexcept -> Pack {
        Pack pack{};
        pack.load(data);
        return pack;
    }

    template<Vector V>
    template<Packet Pack>
    auto CompactVectorView<V>::Iterator::load(size_t count) const noexcept -> Pack {
        assert(0 < count && count < Pack::size() && "[Error]: Invalid size for partial operation");
        Pack pack{};
        pack.load(data, count);
        return pack;
    }

    template<Vector V>
    void CompactVectorView<V>::Iterator::store(const Packet auto pack) noexcept {
        pack.store(data);
    }

    template<Vector V>
    void CompactVectorView<V>::Iterator::store(const Packet auto pack, size_t count) noexcept {
        assert(0 < count && count < pack.size() && "[Error]: Invalid size for partial operation");
        pack.store(data, count);
    }
}
