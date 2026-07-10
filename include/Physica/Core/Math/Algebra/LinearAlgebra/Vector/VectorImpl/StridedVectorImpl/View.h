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

#include <array>
#include "../StridedVector.h"

namespace Physica {
    template<class Derived>
    template<Vector V>
    class StridedVector<Derived>::View final : public std::ranges::view_base {
        using This = View;

        using pointer = decltype(std::declval<V>().data_handle());

        class Iterator;
    private:
        pointer handle;
        size_t stride;
        size_t length;
    public:
        [[gnu::always_inline]] constexpr View(V& vec) noexcept;
        constexpr View(const This&) = default;
        constexpr View(This&&) noexcept = default;
        constexpr ~View() = default;
        /* Operators */
        constexpr This& operator=(const This&) = default;
        constexpr This& operator=(This&&) noexcept = default;
        /* Operations */
        [[nodiscard, gnu::always_inline]] constexpr auto begin(this auto&&) noexcept;
        [[nodiscard, gnu::always_inline]] constexpr auto end(this auto&&) noexcept;
        [[nodiscard, gnu::always_inline]] constexpr size_t size() const noexcept;
    };

    template<class Derived>
    template<Vector V>
    constexpr StridedVector<Derived>::View<V>::View(V& vec) noexcept : handle(vec.data_handle()), stride(vec.getStride()), length(vec.getLength()) {}

    template<class Derived>
    template<Vector V>
    constexpr auto StridedVector<Derived>::View<V>::begin(this auto&& self) noexcept {
        return Iterator(self.handle, self.stride);
    }

    template<class Derived>
    template<Vector V>
    constexpr auto StridedVector<Derived>::View<V>::end(this auto&& self) noexcept {
        return Iterator(self.handle + self.length * self.stride, self.stride);
    }

    template<class Derived>
    template<Vector V>
    constexpr size_t StridedVector<Derived>::View<V>::size() const noexcept {
        return length;
    }
}

namespace Physica {
    template<class Derived>
    template<Vector V>
    class StridedVector<Derived>::View<V>::Iterator {
        using This = Iterator;

        constexpr static size_t StrideAtCompile = Derived::getStrideAtCompile();
        using Stride = std::conditional_t<StrideAtCompile == Dynamic, size_t, Empty>;
    public:
        using iterator_concept = std::random_access_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = V::ScalarType;
        using reference = value_type::RefTy;
        using const_reference = value_type::ConstRefTy;
        using pointer = decltype(std::declval<V>().data_handle());
    private:
        pointer pos;
        Stride stride;
    public:
        constexpr Iterator() = default;
        [[gnu::always_inline]] constexpr Iterator(pointer pos, size_t stride) noexcept;
        constexpr Iterator(const This&) = default;
        constexpr Iterator(This&&) noexcept = default;
        constexpr ~Iterator() = default;
        /* Operators */
        constexpr This& operator=(const This&) = default;
        constexpr This& operator=(This&&) noexcept = default;
        [[gnu::always_inline]] constexpr This& operator++() noexcept;
        [[gnu::always_inline]] constexpr This& operator--() noexcept;
        [[gnu::always_inline]] constexpr This& operator+=(difference_type n) noexcept;
        [[gnu::always_inline]] constexpr This& operator-=(difference_type n) noexcept;
        [[nodiscard, gnu::always_inline]] constexpr This operator++(int) noexcept;
        [[nodiscard, gnu::always_inline]] constexpr This operator--(int) noexcept;
        [[nodiscard, gnu::always_inline]] constexpr decltype(auto) operator*() const noexcept;
        [[nodiscard, gnu::always_inline]] constexpr decltype(auto) operator[](difference_type n) const noexcept;
        [[nodiscard, gnu::always_inline]] constexpr bool operator==(const This& other) const noexcept;
        [[nodiscard, gnu::always_inline]] constexpr auto operator<=>(const This& other) const noexcept;
        [[nodiscard, gnu::always_inline]] constexpr This operator+(difference_type n) const noexcept;
        [[nodiscard, gnu::always_inline]] constexpr This operator-(difference_type n) const noexcept;
        [[nodiscard, gnu::always_inline]] constexpr difference_type operator-(const This& other) const noexcept;
        /* Operations */
        template<int Size>
        [[nodiscard, gnu::always_inline]] SIMD<value_type, Size> load() const noexcept;
        template<int Size>
        [[nodiscard, gnu::always_inline]] SIMD<value_type, Size> load(size_t count) const noexcept;
        [[gnu::always_inline]] void store(Packet auto pack) noexcept;
        [[gnu::always_inline]] void store(Packet auto pack, size_t count) noexcept;
        /* Getters */
        [[nodiscard]] constexpr size_t getStride() const noexcept;
        /* Friends */
        [[gnu::always_inline]] friend constexpr This operator+(difference_type n, const This& ite) noexcept { return ite + n; }
    private:
        template<int Size>
        SIMD<value_type, Size> load_stride2() const noexcept;

        template<int Size>
        SIMD<value_type, Size> load_fallback() const noexcept;
        template<int Size>
        SIMD<value_type, Size> load_fallback(size_t count) const noexcept;
        void store_fallback(Packet auto pack) noexcept;
        void store_fallback(Packet auto pack, size_t count) noexcept;
    };

    template<class Derived>
    template<Vector V>
    constexpr StridedVector<Derived>::View<V>::Iterator::Iterator(pointer pos, size_t stride) noexcept : pos(pos), stride(stride) {}

    template<class Derived>
    template<Vector V>
    constexpr auto StridedVector<Derived>::View<V>::Iterator::operator++() noexcept -> This& {
        pos += getStride();
        return *this;
    }

    template<class Derived>
    template<Vector V>
    constexpr auto StridedVector<Derived>::View<V>::Iterator::operator--() noexcept -> This& {
        pos -= getStride();
        return *this;
    }

    template<class Derived>
    template<Vector V>
    constexpr auto StridedVector<Derived>::View<V>::Iterator::operator+=(difference_type n) noexcept -> This& {
        pos += n * getStride();
        return *this;
    }

    template<class Derived>
    template<Vector V>
    constexpr auto StridedVector<Derived>::View<V>::Iterator::operator-=(difference_type n) noexcept -> This& {
        pos -= n * getStride();
        return *this;
    }

    template<class Derived>
    template<Vector V>
    constexpr auto StridedVector<Derived>::View<V>::Iterator::operator++(int) noexcept -> This {
        return std::exchange(*this, ++(*this));
    }

    template<class Derived>
    template<Vector V>
    constexpr auto StridedVector<Derived>::View<V>::Iterator::operator--(int) noexcept -> This {
        return std::exchange(*this, --(*this));
    }

    template<class Derived>
    template<Vector V>
    constexpr decltype(auto) StridedVector<Derived>::View<V>::Iterator::operator*() const noexcept {
        return *pos;
    }

    template<class Derived>
    template<Vector V>
    constexpr decltype(auto) StridedVector<Derived>::View<V>::Iterator::operator[](difference_type n) const noexcept {
        return *(*this + n);
    }

    template<class Derived>
    template<Vector V>
    constexpr bool StridedVector<Derived>::View<V>::Iterator::operator==(const This& other) const noexcept {
        return pos == other.pos;
    }

    template<class Derived>
    template<Vector V>
    constexpr auto StridedVector<Derived>::View<V>::Iterator::operator<=>(const This& other) const noexcept {
        return pos <=> other.pos;
    }

    template<class Derived>
    template<Vector V>
    constexpr auto StridedVector<Derived>::View<V>::Iterator::operator+(difference_type n) const noexcept -> This {
        return This(pos + n * getStride(), getStride());
    }

    template<class Derived>
    template<Vector V>
    constexpr auto StridedVector<Derived>::View<V>::Iterator::operator-(difference_type n) const noexcept -> This {
        return This(pos - n * getStride(), getStride());
    }

    template<class Derived>
    template<Vector V>
    constexpr auto StridedVector<Derived>::View<V>::Iterator::operator-(const This& other) const noexcept -> difference_type {
        pointer diff = pos - other.pos;
        assert(diff % getStride() == 0 && "[Error]: Corrupted iterator");
        return diff / getStride();
    }

    template<class Derived>
    template<Vector V>
    template<int Size>
    auto StridedVector<Derived>::View<V>::Iterator::load() const noexcept -> SIMD<value_type, Size> {
        if constexpr (StrideAtCompile == 1) {
            SIMD<value_type, Size> pack{};
            pack.load(pos);
            return pack;
        }
        else if constexpr (StrideAtCompile == 2)
            return load_stride2<Size>();
        else
            return load_fallback<Size>();
    }

    template<class Derived>
    template<Vector V>
    template<int Size>
    auto StridedVector<Derived>::View<V>::Iterator::load(size_t count) const noexcept -> SIMD<value_type, Size> {
        assert(0 < count && count < Size && "[Error]: Invalid size for partial operation");
        if constexpr (StrideAtCompile == 1) {
            SIMD<value_type, Size> pack{};
            pack.load(pos, count);
            return pack;
        }
        else
            return load_fallback<Size>(count);
    }

    template<class Derived>
    template<Vector V>
    void StridedVector<Derived>::View<V>::Iterator::store(const Packet auto pack) noexcept {
        if constexpr (StrideAtCompile == 1)
            pack.store(pos);
        else
            store_fallback(pack);
    }

    template<class Derived>
    template<Vector V>
    void StridedVector<Derived>::View<V>::Iterator::store(const Packet auto pack, size_t count) noexcept {
        assert(0 < count && count < pack.size() && "[Error]: Invalid size for partial operation");
        if constexpr (StrideAtCompile == 1)
            pack.store(pos, count);
        else
            store_fallback(pack, count);
    }

    template<class Derived>
    template<Vector V>
    template<int Size>
    auto StridedVector<Derived>::View<V>::Iterator::load_stride2() const noexcept -> SIMD<value_type, Size> {
        using Pack = SIMD<value_type, Size>;
        constexpr int MaxSize = BestPacket<value_type, Dynamic>::Size;
        if constexpr (2 * Size <= MaxSize) {
            SIMD<value_type, 2 * Size> pack2;
            pack2.load(pos);
            pack2.gatherRealImag();
            return pack2.getLow();
        }
        else {
            Pack low, high;
            low.load(pos);
            high.load(pos + Size * getStride());
            return Pack(low.gatherRealImag().getLow(), high.gatherRealImag().getLow());
        }
    }

    template<class Derived>
    template<Vector V>
    template<int Size>
    auto StridedVector<Derived>::View<V>::Iterator::load_fallback() const noexcept -> SIMD<value_type, Size> {
        Array<value_type, Size> buffer{};
        for (size_t i = 0; i < Size; ++i)
            buffer[i] = (*this)[i];
        SIMD<value_type, Size> pack{};
        pack.load(buffer.data());
        return pack;
    }

    template<class Derived>
    template<Vector V>
    template<int Size>
    auto StridedVector<Derived>::View<V>::Iterator::load_fallback(size_t count) const noexcept -> SIMD<value_type, Size> {
        Array<value_type, Size> buffer{};
        for (size_t i = 0; i < count; ++i)
            buffer[i] = (*this)[i];
        SIMD<value_type, Size> pack{};
        pack.load(buffer.data(), count);
        return pack;
    }

    template<class Derived>
    template<Vector V>
    void StridedVector<Derived>::View<V>::Iterator::store_fallback(const Packet auto pack) noexcept {
        Array<value_type, pack.size()> buffer{};
        pack.store(buffer.data());
        for (size_t i = 0; i < pack.size(); ++i)
            (*this)[i] = buffer[i];
    }

    template<class Derived>
    template<Vector V>
    void StridedVector<Derived>::View<V>::Iterator::store_fallback(const Packet auto pack, size_t count) noexcept {
        assert(0 < count && count < pack.size() && "[Error]: Invalid size for partial operation");
        Array<value_type, pack.size()> buffer{};
        pack.store(buffer.data(), count);
        for (size_t i = 0; i < count; ++i)
            (*this)[i] = buffer[i];
    }

    template<class Derived>
    template<Vector V>
    constexpr size_t StridedVector<Derived>::View<V>::Iterator::getStride() const noexcept {
        if constexpr (StrideAtCompile == Dynamic)
            return stride;
        else
            return StrideAtCompile;
    }
}
