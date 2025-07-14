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

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <limits>
#include "Physica/Macro.h"

namespace Physica {
    /**
     * Default allocator for \class Array, which provides custom interface reallocate()
     * 
     * Designed to meet C++17 standard(WIP)
     */
    template<class T, size_t Align = Dynamic>
    class HostAllocator {
        using This = HostAllocator<T, Align>;

        constexpr static bool OverAlign = Align > alignof(std::max_align_t);
    public:
        using value_type = T;
        using size_type = std::size_t;
        using difference_type = std::ptrdiff_t;
        using propagate_on_container_move_assignment = std::true_type;
        template<class U>
        using rebind_alloc = HostAllocator<U, Align>;
    public:
        HostAllocator() noexcept = default;
        HostAllocator(const This&) noexcept = default;
        HostAllocator(This&&) noexcept = delete;
        ~HostAllocator() = default;
        /* Operators */
        This& operator=(const This&) noexcept = default;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T* allocate(size_t n) noexcept;
        void deallocate(T* p, size_t n) noexcept;
        [[nodiscard]] T* reallocate(T* p, size_t new_size, size_t old_size) noexcept;
        void construct(T* p, auto&&... args);
        void destroy(T* p) noexcept;
    };

    template<class T, size_t Align>
    [[nodiscard]] T* HostAllocator<T, Align>::allocate(size_t n) noexcept {
        assert(n > 0 && "[Error]: Allocate nothing");
        size_t size = n * sizeof(T);
        void* p;
        if constexpr (OverAlign) {
            size = (size + Align - 1) / Align * Align;
        #ifdef _MSC_VER
            p = _aligned_malloc(Align, size);
        #else
            p = std::aligned_alloc(Align, size);
        #endif
        }
        else
            p = std::malloc(size);
        return reinterpret_cast<T*>(p);
    }

    template<class T, size_t Align>
    void HostAllocator<T, Align>::deallocate(T* p, [[maybe_unused]] size_t n) noexcept {
        if constexpr (OverAlign) {
        #ifdef _MSC_VER
            _aligned_free(p);
        #else
            std::free(p);
        #endif
        }
        else
            std::free(p);
    }

    template<class T, size_t Align>
    T* HostAllocator<T, Align>::reallocate(T* p, size_t new_size, [[maybe_unused]] size_t old_size) noexcept {
        assert(new_size > 0 && "[Error]: Allocate nothing");
        T* new_p = allocate(new_size);
        if (p == nullptr)
            return new_p;

        const size_t size = std::min(new_size, old_size);
        if constexpr (std::is_trivially_copyable<T>::value)
            memcpy(new_p, p, size * sizeof(T));
        else {
            for (size_t i = 0; i < size; ++i)
                construct(new_p + i, std::move(p[i]));
        }
        deallocate(p, old_size);
        return new_p;
    }

    template<class T, size_t Align>
    void HostAllocator<T, Align>::construct(T* p, auto&&... args) {
        ::new (static_cast<void*>(p)) T(std::forward<decltype(args)>(args)...);
    }

    template<class T, size_t Align>
    void HostAllocator<T, Align>::destroy(T* p) noexcept {
        p->~T();
    }
}

namespace std {
    template<class T, size_t Align_>
    struct allocator_traits<Physica::HostAllocator<T, Align_>> {
    public:
        using allocator_type = Physica::HostAllocator<T, Align_>;
        using value_type = T;
        using pointer = T*;
        using const_pointer = const T*;
        using void_pointer = void*;
        using const_void_pointer = const void*;
        using lvalue_reference = T&;
        using const_lvalue_reference = const T&;
        using rvalue_reference = T&&;
        using size_type = allocator_type::size_type;
        using difference_type = allocator_type::difference_type;
        using propagate_on_container_copy_assignment = std::false_type;
        using propagate_on_container_move_assignment = std::false_type;
        using propagate_on_container_swap = std::false_type;
        using is_always_equal = std::is_empty<allocator_type>::type;
        template<class U>
        using rebind_alloc = Physica::HostAllocator<U>;
        template<class U>
        using rebind_traits = std::allocator_traits<rebind_alloc<U>>;

        constexpr static size_t Align = Align_;
        constexpr static bool isPageLocked = false;
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

        static void construct(allocator_type& a, T* p, auto&&... args) {
            a.construct(p, std::forward<decltype(args)>(args)...);
        }

        static void destroy(allocator_type& a, T* p) {
            a.destroy(p);
        }

        static constexpr size_type max_size([[maybe_unused]] const allocator_type& a) noexcept {
            return std::numeric_limits<size_type>::max() / sizeof(value_type);
        }

        static allocator_type select_on_container_copy_construction(const allocator_type& a) {
            allocator_type result = a;
            return result;
        }
    };
}
