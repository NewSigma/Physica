/*
 * Copyright 2021-2024 Weibo He.
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
#include <memory>
#include "Physica/Core/Parallel/CUDAContext.cuh"
#include "Physica/Core/Utils/CUDA/device_obj.cuh"

namespace Physica {
    namespace Internal {
        template<class T, bool IsClass>
        struct DeviceAllocatorValueType {
            using Type = T::device_obj_type;
        };

        template<class T>
        struct DeviceAllocatorValueType<T, false> {
            using Type = T;
        };
    }

    template<class T>
    class DeviceAllocator {
    public:
        using value_type = Internal::DeviceAllocatorValueType<T, std::is_class<T>::value>::Type;
        using pointer = value_type*;
        using size_type = size_t;
        using difference_type = std::ptrdiff_t;
        using propagate_on_container_move_assignment = std::true_type;
        template<class U>
        using rebind_alloc = DeviceAllocator<U>;
    public:
        DeviceAllocator() noexcept = default;
        DeviceAllocator(const DeviceAllocator&) noexcept = default;
        DeviceAllocator(DeviceAllocator&&) noexcept = default;
        ~DeviceAllocator() = default;
        /* Operators */
        DeviceAllocator& operator=(const DeviceAllocator&) noexcept = default;
        DeviceAllocator& operator=(DeviceAllocator&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __host__ __device__ pointer allocate(size_t n) noexcept;
        __host__ __device__ void deallocate(pointer p, size_t n) noexcept;
        __host__ __device__ void construct(pointer p, auto&&... args);
        __host__ __device__ void destroy(pointer p) noexcept;
    };

    template<class T>
    __host__ __device__ DeviceAllocator<T>::pointer DeviceAllocator<T>::allocate(size_t n) noexcept {
        pointer p;
        if constexpr (IsDevice())
            p = reinterpret_cast<pointer>(malloc(n * sizeof(value_type)));
        else {
            if constexpr (CUDADevAttr::MemoryPoolsSupported)
                cudaMallocAsync(&p, n * sizeof(value_type), CUDAContext::getInstance());
            else
                cudaMalloc(&p, n * sizeof(value_type));
        }
        return p;
    }

    template<class T>
    __host__ __device__ void DeviceAllocator<T>::deallocate(pointer p, [[maybe_unused]] size_t n) noexcept {
        if constexpr (IsDevice())
            std::free(p);
        else {
            if (p != nullptr) { // No unnecessary cuda api call to make profiler output cleaner
                if constexpr (CUDADevAttr::MemoryPoolsSupported)
                    cudaFreeAsync(p, CUDAContext::getInstance());
                else
                    cudaFree(p);
            }
        }
    }

    template<class T>
    __host__ __device__ void DeviceAllocator<T>::construct(pointer p, auto&&... args) {
        if constexpr (IsDevice())
            ::new (static_cast<void*>(p)) value_type(std::forward<decltype(args)>(args)...);
        else {
            value_type temp(std::forward<decltype(args)>(args)...);
            check(cudaMemcpyAsync(p, &temp, sizeof(value_type), cudaMemcpyHostToDevice, CUDAContext::getInstance()));
            if constexpr (!std::is_trivially_copyable<value_type>::value)
                temp.release(); //Ownership has been given to device
        }
    }

    template<class T>
    __host__ __device__ void DeviceAllocator<T>::destroy(pointer p) noexcept {
        if constexpr (!std::is_trivially_copyable<value_type>::value) {
            if constexpr (IsDevice())
                p->~value_type();
            else {
                value_type temp;
                cudaMemcpyAsync(&temp, p, sizeof(value_type), cudaMemcpyDeviceToHost, CUDAContext::getInstance());
            }
        }
    }
}

namespace std {
    template<class T>
    struct allocator_traits<Physica::DeviceAllocator<T>> {
    public:
        using allocator_type = Physica::DeviceAllocator<T>;
        using value_type = allocator_type::value_type;
        using pointer = value_type*;
        using const_pointer = const value_type*;
        using void_pointer = void*;
        using const_void_pointer = const void*;
        using lvalue_reference = value_type&;
        using const_lvalue_reference = const value_type&;
        using rvalue_reference = value_type&&;
        using size_type = size_t;
        using difference_type = std::ptrdiff_t;
        using propagate_on_container_copy_assignment = std::false_type;
        using propagate_on_container_move_assignment = std::false_type;
        using propagate_on_container_swap = std::false_type;
        using is_always_equal = std::is_empty<allocator_type>::type;
        template<class U>
        using rebind_alloc = Physica::DeviceAllocator<U>;
        template<class U>
        using rebind_traits = std::allocator_traits<rebind_alloc<U>>;

        constexpr static size_t Align = 256; // CUDA default
        constexpr static bool isPageLocked = false;
    public:
        [[nodiscard]] __host__ __device__ static pointer allocate(allocator_type& a, size_type n) {
            return a.allocate(n);
        }

        __host__ __device__ static void deallocate(allocator_type& a, pointer p, size_type n) noexcept {
            a.deallocate(p, n);
        }

        __host__ __device__ static void construct(allocator_type& a, pointer p, auto&&... args) {
            a.construct(p, std::forward<decltype(args)>(args)...);
        }

        __host__ __device__ static void destroy(allocator_type& a, pointer p) {
            a.destroy(p);
        }

         __host__ __device__ static constexpr size_type max_size([[maybe_unused]] const allocator_type& a) noexcept {
            return std::numeric_limits<size_type>::max() / sizeof(value_type);
        }

         __host__ __device__ static allocator_type select_on_container_copy_construction(const allocator_type& a) {
            allocator_type result = a;
            return result;
        }
    };
}
