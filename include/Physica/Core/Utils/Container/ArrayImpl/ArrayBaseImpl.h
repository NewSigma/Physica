/*
 * Copyright 2020-2026 Weibo He.
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
#include <cstdio>
#include <cstring>
#include "ArrayBase.h"
#include "Physica/Core/IO/HDF5/HDF5.h"

namespace Physica {
    template<class Derived, class Allocator>
    __host__ __device__ constexpr auto& ArrayBase<Derived, Allocator>::operator[](this auto&& self, size_t index) noexcept {
        assert(index < self.getLength() && "[Error]: Index overflow");
        return self.data()[index];
    }

    template<class Derived, class Allocator>
    __host__ __device__ bool ArrayBase<Derived, Allocator>::operator==(this const This& lhs, const This& rhs) {
        if (lhs.getLength() != rhs.getLength())
            return false;
        if (lhs.getCapacity() != rhs.getCapacity())
            return false;
        for (size_t i = 0; i < lhs.getLength(); ++i)
            if (lhs[i] != rhs[i])
                return false;
        return true;
    }

    template<class Derived, class Allocator>
    __host__ __device__ bool ArrayBase<Derived, Allocator>::operator!=(this const This& lhs, const This& rhs) {
        return !(lhs == rhs);
    }

    template<class Derived, class Allocator>
    __host__ __device__ constexpr auto ArrayBase<Derived, Allocator>::begin(this auto&& self) noexcept {
        using Self = decltype(self);
        using Container = std::conditional<std::is_const_v<std::remove_reference_t<Self>>, const Derived, Derived>::type;
        return Iterator<Container>(self.data());
    }

    template<class Derived, class Allocator>
    __host__ __device__ constexpr auto ArrayBase<Derived, Allocator>::end(this auto&& self) noexcept {
        using Self = decltype(self);
        using Container = std::conditional<std::is_const_v<std::remove_reference_t<Self>>, const Derived, Derived>::type;
        return Iterator<Container>(self.data() + self.getLength());
    }

    template<class Derived, class Allocator>
    template<class R>
    auto& ArrayBase<Derived, Allocator>::select(this auto&& self) noexcept {
        assert(!self.empty() && "[Error]: Cannot select from empty array");
        return self[R::select(self.getLength())];
    }

    template<class Derived, class Allocator>
    void ArrayBase<Derived, Allocator>::read(this auto& self, const auto& loc, const char* name) {
        if constexpr (std::is_trivially_copyable_v<value_type>) {
            const auto dataset = loc.template openDataSet<1>(name);
            const size_t length = dataset.getSize(0);
            const auto space = H5DataSpace<1>(length);
            self.resize(length);
            dataset.read(self.data(), H5Type::get<value_type>(), space, space);
        }
        else {
            auto group = loc.openGroup(name);
            std::array<char, 32> buffer{}; // 32 is enough for uint64_t
            for (size_t i = 0; i < self.getLength(); ++i) {
                [[maybe_unused]] const auto count = std::format_to_n(buffer.data(), buffer.size(), "%zu", i).size;
                assert(count <= buffer.size());
                self[i].read(group, buffer.data());
            }
        }
    }

    template<class Derived, class Allocator>
    void ArrayBase<Derived, Allocator>::write(this const auto& self, auto& loc, const char* name) {
        if constexpr (std::is_trivially_copyable_v<value_type>) {
            const size_t length = self.getLength();
            const auto space = H5DataSpace<1>(length);
            H5Dataset<1> dataset;
            if (loc.exists(name))
                dataset = loc.template openDataSet<1>(name);
            else
                dataset = loc.template createDataSet<1>(name, H5Type::get<value_type>(), space);
            dataset.write(self.data(), H5Type::get<value_type>(), space, space);
        }
        else {
            auto group = loc.openGroup(name);
            std::array<char, 32> buffer{}; // 32 is enough for uint64_t
            for (size_t i = 0; i < self.getLength(); ++i) {
                [[maybe_unused]] const auto count = std::format_to_n(buffer.data(), buffer.size(), "%zu", i).size;
                assert(count <= buffer.size());
                self[i].write(group, buffer.data());
            }
        }
    }

    template<class Derived, class Allocator>
    __host__ __device__ auto* ArrayBase<Derived, Allocator>::data(this auto&& self) noexcept {
        return self.getDerived().data();
    }

    template<class Derived, class Allocator>
    __host__ __device__ auto* ArrayBase<Derived, Allocator>::data_ptr(this auto&& self, size_t index) noexcept {
        assert(index < self.getLength());
        return self.data() + index;
    }

    template<class Derived, class Allocator>
    __host__ __device__ auto& ArrayBase<Derived, Allocator>::front(this auto&& self) noexcept {
        assert(!self.empty());
        return self[0];
    }

    template<class Derived, class Allocator>
    __host__ __device__ auto& ArrayBase<Derived, Allocator>::back(this auto&& self) noexcept {
        assert(!self.empty());
        return self[self.getLength() - 1];
    }

    template<class Derived, class Allocator>
    template<class... Args>
    __host__ __device__ consteval bool ArrayBase<Derived, Allocator>::isTrivialDefaultConstruct() noexcept {
        return (sizeof...(Args) == 0) && std::is_trivially_default_constructible<value_type>::value;
    }

    template<class Derived, class Allocator>
    template<class Container>
    __host__ __device__ constexpr ArrayBase<Derived, Allocator>::Iterator<Container>::Iterator(pointer p) noexcept : p(p) {}

    template<class Derived, class Allocator>
    template<class Container>
    __host__ __device__ constexpr auto ArrayBase<Derived, Allocator>::Iterator<Container>::operator++() noexcept -> This& {
        ++p;
        return *this;
    }

    template<class Derived, class Allocator>
    template<class Container>
    __host__ __device__ constexpr auto ArrayBase<Derived, Allocator>::Iterator<Container>::operator--() noexcept -> This& {
        --p;
        return *this;
    }

    template<class Derived, class Allocator>
    template<class Container>
    __host__ __device__ constexpr auto ArrayBase<Derived, Allocator>::Iterator<Container>::operator+=(difference_type n) noexcept -> This& {
        p += n;
        return *this;
    }

    template<class Derived, class Allocator>
    template<class Container>
    __host__ __device__ constexpr auto ArrayBase<Derived, Allocator>::Iterator<Container>::operator-=(difference_type n) noexcept -> This& {
        p -= n;
        return *this;
    }

    template<class Derived, class Allocator>
    template<class Container>
    __host__ __device__ constexpr auto ArrayBase<Derived, Allocator>::Iterator<Container>::operator++(int) noexcept -> This {
        return This(p++);
    }

    template<class Derived, class Allocator>
    template<class Container>
    __host__ __device__ constexpr auto ArrayBase<Derived, Allocator>::Iterator<Container>::operator--(int) noexcept -> This {
        return This(p--);
    }

    template<class Derived, class Allocator>
    template<class Container>
    __host__ __device__ constexpr auto ArrayBase<Derived, Allocator>::Iterator<Container>::operator*() const noexcept -> reference {
        return *p;
    }

    template<class Derived, class Allocator>
    template<class Container>
    __host__ __device__ constexpr auto ArrayBase<Derived, Allocator>::Iterator<Container>::operator->() const noexcept -> pointer {
        return p;
    }

    template<class Derived, class Allocator>
    template<class Container>
    __host__ __device__ constexpr auto ArrayBase<Derived, Allocator>::Iterator<Container>::operator[](difference_type n) const noexcept -> reference {
        return *operator+(n);
    }

    template<class Derived, class Allocator>
    template<class Container>
    __host__ __device__ constexpr bool ArrayBase<Derived, Allocator>::Iterator<Container>::operator==(const This& other) const noexcept {
        return p == other.p;
    }

    template<class Derived, class Allocator>
    template<class Container>
    __host__ __device__ constexpr auto ArrayBase<Derived, Allocator>::Iterator<Container>::operator<=>(const This& other) const noexcept {
        return p <=> other.p;
    }

    template<class Derived, class Allocator>
    template<class Container>
    __host__ __device__ constexpr auto ArrayBase<Derived, Allocator>::Iterator<Container>::operator+(difference_type n) const noexcept -> This {
        return This(p + n);
    }

    template<class Derived, class Allocator>
    template<class Container>
    __host__ __device__ constexpr auto ArrayBase<Derived, Allocator>::Iterator<Container>::operator-(difference_type n) const noexcept -> This {
        return This(p - n);
    }

    template<class Derived, class Allocator>
    template<class Container>
    __host__ __device__ constexpr auto ArrayBase<Derived, Allocator>::Iterator<Container>::operator-(const This& other) const noexcept -> difference_type {
        return p - other.p;
    }
}
