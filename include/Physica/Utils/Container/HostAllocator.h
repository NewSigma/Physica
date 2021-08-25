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

#include <cstdlib>
#include <new>

namespace Physica::Utils {
    /**
     * Default allocator for \class Array, which provides custom interface reallocate()
     * 
     * Designed to meet C++17 standard(WIP)
     */
    template<class T>
    class HostAllocator {
    public:
        using value_type = T;
        using size_type = std::size_t;
        using difference_type = std::ptrdiff_t;
        using propagate_on_container_move_assignment = std::true_type;
    public:
        HostAllocator() noexcept = default;
        HostAllocator(const HostAllocator&) noexcept = default;
        HostAllocator(HostAllocator&&) noexcept = delete;
        ~HostAllocator() = default;
        /* Operators */
        HostAllocator& operator=(const HostAllocator&) noexcept = default;
        HostAllocator& operator=(HostAllocator&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] static T* allocate(size_t n);
        static void deallocate(T* p, size_t n);
        [[nodiscard]] static T* reallocate(T* p, size_t n);
    };

    template<class T>
    T* HostAllocator<T>::allocate(size_t n) {
        auto* p = reinterpret_cast<T*>(malloc(n * sizeof(T)));
        if (!p)
            throw std::bad_alloc();
        return p;
    }

    template<class T>
    void HostAllocator<T>::deallocate(T* p, [[maybe_unused]] size_t n) {
        free(p);
    }

    template<class T>
    T* HostAllocator<T>::reallocate(T* p, size_t n) {
        return reinterpret_cast<T*>(realloc(p, n * sizeof(T)));
    }
}
