/*
 * Copyright 2023 WeiBo He.
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

#include "Physica/Utils/CUDA/DebugUtil.cuh"
#include "HostAllocator.h"

namespace Physica::Utils {
    template<class T> class PageLockedAllocator;
}

namespace std {
    template<class T>
    struct allocator_traits<Physica::Utils::PageLockedAllocator<T>> {
    public:
        using allocator_type = Physica::Utils::PageLockedAllocator<T>;
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
        using is_always_equal = typename std::is_empty<allocator_type>::type;
        template<class U>
        using rebind_alloc = Physica::Utils::PageLockedAllocator<U>;
        template<class U>
        using rebind_traits = std::allocator_traits<rebind_alloc<U>>;
    public:
        [[nodiscard]] __host__ __device__ static pointer allocate(allocator_type& a, size_type n) {
            return a.allocate(n);
        }

        __host__ __device__ static void deallocate(allocator_type& a, pointer p, size_type n) noexcept {
            a.deallocate(p, n);
        }

        template<class... Args>
        __host__ __device__ static void construct(allocator_type& a, pointer p, Args&&... args) {
            a.construct(p, std::forward<Args>(args)...);
        }

        __host__ __device__ static void destroy(allocator_type& a, pointer p) {
            a.destroy(p);
        }

         __host__ __device__ static constexpr size_type max_size(const allocator_type& a) noexcept {
            return std::numeric_limits<size_type>::max() / sizeof(value_type);
        }

         __host__ __device__ static allocator_type select_on_container_copy_construction(const allocator_type& a) {
            allocator_type result = a;
            return result;
        }
    };
}

namespace Physica::Utils {
    template<class T>
    class PageLockedAllocator {
    public:
        using value_type = typename std::allocator_traits<PageLockedAllocator>::value_type;
        using pointer = typename std::allocator_traits<PageLockedAllocator>::pointer;
        using size_type = typename std::allocator_traits<PageLockedAllocator>::size_type;
        using difference_type = typename std::allocator_traits<PageLockedAllocator>::difference_type;
        using propagate_on_container_move_assignment = std::true_type;
    public:
        PageLockedAllocator() noexcept = default;
        PageLockedAllocator(const PageLockedAllocator&) noexcept = default;
        PageLockedAllocator(PageLockedAllocator&&) noexcept = default;
        ~PageLockedAllocator() = default;
        /* Operators */
        PageLockedAllocator& operator=(const PageLockedAllocator&) noexcept = default;
        PageLockedAllocator& operator=(PageLockedAllocator&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] static pointer allocate(size_t n);
        inline static void deallocate(pointer p, size_t n) noexcept;
        template<class... Args>
        inline void construct(pointer p, Args&&... args);
        inline void destroy(pointer p);
    };

    template<class T>
    typename PageLockedAllocator<T>::pointer PageLockedAllocator<T>::allocate(size_t n) {
        pointer p;
        cudaCheck(cudaMallocHost(&p, n * sizeof(value_type)));
        return p;
    }

    template<class T>
    inline void PageLockedAllocator<T>::deallocate(pointer p, [[maybe_unused]] size_t n) noexcept {
        cudaFreeHost(p);
    }

    template<class T>
    template<class... Args>
    inline void PageLockedAllocator<T>::construct(pointer p, Args&&... args) {
        HostAllocator<T> alloc{};
        alloc.construct(p, std::forward<Args>(args)...);
    }

    template<class T>
    inline void PageLockedAllocator<T>::destroy(pointer p) {
        HostAllocator<T> alloc{};
        alloc.destroy(p);
    }
}
