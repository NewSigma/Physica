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

#include "HostAllocator.h"

namespace Physica {
    template<class T>
    class PageLockedAllocator {
        using This = PageLockedAllocator<T>;
    public:
        using value_type = T;
        using pointer = T*;
        using size_type = size_t;
        using difference_type = std::ptrdiff_t;
        using propagate_on_container_move_assignment = std::true_type;
        template<class U>
        using rebind_alloc = PageLockedAllocator<U>;
    public:
        PageLockedAllocator() noexcept = default;
        PageLockedAllocator(const This&) noexcept = default;
        PageLockedAllocator(This&&) noexcept = default;
        ~PageLockedAllocator() = default;
        /* Operators */
        This& operator=(const This&) noexcept = default;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard, gnu::returns_nonnull]] static pointer allocate(size_t n);
        static void deallocate(pointer p, size_t n) noexcept;
        [[nodiscard, gnu::returns_nonnull]] T* reallocate(T* p, size_t new_size, size_t old_size);
        void construct(pointer p, auto&&... args);
        void destroy(pointer p);
    };

    template<class T>
    auto PageLockedAllocator<T>::allocate(size_t n) -> pointer {
        assert(n > 0 && "[Error]: Allocate nothing");
        pointer p{};
        check(cudaMallocHost(&p, n * sizeof(value_type)));
        return p;
    }

    template<class T>
    void PageLockedAllocator<T>::deallocate(pointer p, [[maybe_unused]] size_t n) noexcept {
        if (p != nullptr)
            cudaFreeHost(p);
    }

    template<class T>
    T* PageLockedAllocator<T>::reallocate(T* p, size_t new_size, [[maybe_unused]] size_t old_size) {
        assert(new_size > 0 && "[Error]: Allocate nothing");
        T* new_p = allocate(new_size);
        size_t size = std::min(new_size, old_size);
        if constexpr (std::is_trivially_copyable<T>::value)
            memcpy(new_p, p, size * sizeof(T));
        else {
            for (size_t i = 0; i < size; ++i)
                construct(new_p + i, std::move(p[i]));
        }
        deallocate(p, old_size);
        return new_p;
    }

    template<class T>
    void PageLockedAllocator<T>::construct(pointer p, auto&&... args) {
        HostAllocator<T> alloc{};
        alloc.construct(p, std::forward<decltype(args)>(args)...);
    }

    template<class T>
    void PageLockedAllocator<T>::destroy(pointer p) {
        HostAllocator<T> alloc{};
        alloc.destroy(p);
    }
}

namespace std {
    template<class T>
    struct allocator_traits<Physica::PageLockedAllocator<T>> {
    public:
        using allocator_type = Physica::PageLockedAllocator<T>;
        using value_type = T;
        using pointer = T*;
        using const_pointer = const T*;
        using void_pointer = void*;
        using const_void_pointer = const void*;
        using lvalue_reference = T&;
        using const_lvalue_reference = const T&;
        using rvalue_reference = T&&;
        using size_type = size_t;
        using difference_type = std::ptrdiff_t;
        using propagate_on_container_copy_assignment = std::false_type;
        using propagate_on_container_move_assignment = std::false_type;
        using propagate_on_container_swap = std::false_type;
        using is_always_equal = std::is_empty<allocator_type>::type;
        template<class U>
        using rebind_alloc = Physica::PageLockedAllocator<U>;
        template<class U>
        using rebind_traits = std::allocator_traits<rebind_alloc<U>>;

        constexpr static size_t Align = 256; // CUDA default
        constexpr static bool isPageLocked = true;
    public:
        [[nodiscard]] static pointer allocate(allocator_type& a, size_type n) {
            return a.allocate(n);
        }

        static void deallocate(allocator_type& a, pointer p, size_type n) noexcept {
            a.deallocate(p, n);
        }

        [[nodiscard]] static pointer reallocate(allocator_type& a, pointer p, size_type n) {
            return a.reallocate(p, n);
        }

        static void construct(allocator_type& a, pointer p, auto&&... args) {
            a.construct(p, std::forward<decltype(args)>(args)...);
        }

        static void destroy(allocator_type& a, pointer p) {
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
