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
#include <memory>

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
        [[nodiscard]] T* allocate(size_t n);
        void deallocate(T* p, size_t n) noexcept;
        [[nodiscard]] T* reallocate(T* p, size_t new_size, size_t old_size);
        template<class... Args>
        void construct(T* p, Args&&... args);
        void destroy(T* p);
    };

    template<class T>
    [[nodiscard]] T* HostAllocator<T>::allocate(size_t n) {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Walloc-size-larger-than="
        auto* p = reinterpret_cast<T*>(malloc(n * sizeof(T)));
    #pragma GCC diagnostic pop
        if (!p)
            throw std::bad_alloc();
        return p;
    }

    template<class T>
    void HostAllocator<T>::deallocate(T* p, [[maybe_unused]] size_t n) noexcept {
        free(p);
    }

    template<class T>
    [[nodiscard]] T* HostAllocator<T>::reallocate(T* p, size_t new_size, [[maybe_unused]] size_t old_size) {
        if constexpr (std::is_trivially_copyable<T>::value)
            return reinterpret_cast<T*>(realloc(p, new_size * sizeof(T)));
        else {
            T* new_p = allocate(new_size);
            for (size_t i = 0; i < std::min(new_size, old_size); ++i)
                construct(new_p + i, std::move(p[i]));
            deallocate(p, old_size);
            return new_p;
        }
    }

    template<class T>
    template<class... Args>
    void HostAllocator<T>::construct(T* p, Args&&... args) {
        ::new (static_cast<void*>(p)) T(std::forward<Args>(args)...);
    }

    template<class T>
    void HostAllocator<T>::destroy(T* p) {
        p->~T();
    }
}

namespace std {
    template<class T>
    struct allocator_traits<Physica::Utils::HostAllocator<T>> {
    public:
        using allocator_type = Physica::Utils::HostAllocator<T>;
        using value_type = T;
        using pointer = T*;
        using const_pointer = const T*;
        using void_pointer = void*;
        using const_void_pointer = const void*;
        using lvalue_reference = T&;
        using const_lvalue_reference = const T&;
        using rvalue_reference = T&&;
        using size_type = typename allocator_type::size_type;
        using difference_type = typename allocator_type::difference_type;
        using propagate_on_container_copy_assignment = std::false_type;
        using propagate_on_container_move_assignment = std::false_type;
        using propagate_on_container_swap = std::false_type;
        using is_always_equal = typename std::is_empty<allocator_type>::type;
        template<class U>
        using rebind_alloc = Physica::Utils::HostAllocator<U>;
        template<class U>
        using rebind_traits = std::allocator_traits<rebind_alloc<U>>;
    public:
        [[nodiscard]] static pointer allocate(allocator_type& a, size_type n) {
            return a.allocate(n);
        }

        static void deallocate(allocator_type& a, pointer p, size_type n) {
            a.deallocate(p, n);
        }

        [[nodiscard]] static pointer reallocate(allocator_type& a, pointer p, size_type n) {
            return a.reallocate(p, n);
        }

        template<class... Args>
        static void construct(allocator_type& a, T* p, Args&&... args) {
            a.construct(p, std::forward<Args>(args)...);
        }

        static void destroy(allocator_type& a, T* p) {
            a.destroy(p);
        }

        static constexpr size_type max_size(const allocator_type& a) noexcept {
            return std::numeric_limits<size_type>::max() / sizeof(value_type);
        }

        static allocator_type select_on_container_copy_construction(const allocator_type& a) {
            allocator_type result = a;
            return result;
        }
    };
}
