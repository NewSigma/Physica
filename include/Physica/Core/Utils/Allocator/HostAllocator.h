/*
 * Copyright 2021-2026 Weibo He.
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
     * Default allocator in Physica
     */
    template<class T, size_t Align = Dynamic>
    class HostAllocator {
        using This = HostAllocator<T, Align>;

        constexpr static bool OverAlign = Align > alignof(std::max_align_t);
    public:
        using value_type = T;
    public:
        HostAllocator() = default;
        HostAllocator(const This&) = default;
        HostAllocator(This&&) noexcept = default;
        ~HostAllocator() = default;
        /* Operators */
        This& operator=(const This&) noexcept = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        [[nodiscard, gnu::returns_nonnull]] T* allocate(size_t n) noexcept;
        void deallocate(T* p, size_t n) noexcept;
        [[nodiscard, gnu::returns_nonnull]] T* reallocate(T* p, size_t new_size, size_t old_size) noexcept;
    private:
        T* reallocate_mimalloc(T* p, size_t new_size, size_t old_size) noexcept;
    };

    template<class T, size_t Align>
    T* HostAllocator<T, Align>::allocate(size_t n) noexcept {
        size_t size = n * sizeof(T);
        void* p{};
        if constexpr (OverAlign)
            p = ::operator new(size, std::align_val_t(Align), std::nothrow);
        else
            p = ::operator new(size, std::nothrow);
        assert(p != nullptr); // null return value is rare in reality
        return reinterpret_cast<T*>(p);
    }

    template<class T, size_t Align>
    void HostAllocator<T, Align>::deallocate(T* p, size_t) noexcept {
        if constexpr (OverAlign)
            ::operator delete(p, std::align_val_t(Align));
        else
            ::operator delete(p);
    }
    /**
     * Reference:
     * [1] N3322; https://www.open-std.org/jtc1/sc22/wg14/www/docs/n3322.pdf
     */
    template<class T, size_t Align>
    T* HostAllocator<T, Align>::reallocate(T* p, size_t new_size, size_t old_size) noexcept {
        if constexpr (HasMimalloc())
            return reallocate_mimalloc(p, new_size, old_size);
        else {
            assert(new_size > 0 && "[Error]: Reject bad pattern");
            assert(p != nullptr || old_size == 0); // According to [1], the behavior is well defined now
            T* new_p = allocate(new_size);
            memcpy((void*)new_p, p, std::min(new_size, old_size) * sizeof(T));
            deallocate(p, old_size);
            return new_p;
        }
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
        using size_type = std::size_t;
        using difference_type = std::ptrdiff_t;
        using propagate_on_container_copy_assignment = std::false_type;
        using propagate_on_container_move_assignment = std::false_type;
        using propagate_on_container_swap = std::false_type;
        using is_always_equal = std::is_empty<allocator_type>::type;
        template<class U>
        using rebind_alloc = Physica::HostAllocator<U, Align_>;
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

        static void construct(allocator_type&, T* p, auto&&... args) noexcept(std::is_nothrow_constructible<T, decltype(args)...>::value) {
            ::new (static_cast<void*>(p)) T(std::forward<decltype(args)>(args)...);
        }

        static void destroy(allocator_type&, T* p) {
            p->~T();
        }

        constexpr static size_type max_size(const allocator_type&) noexcept {
            return std::numeric_limits<size_type>::max() / sizeof(value_type);
        }

        static allocator_type select_on_container_copy_construction(const allocator_type& a) {
            allocator_type result = a;
            return result;
        }
    };
}

#ifdef PHYSICA_MIMALLOC
    #include "HostAllocator_Mimalloc.h"
#endif
