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

#include <cstdlib>
#include <memory>
#include "Physica/Core/Parallel/CUDAContext.cuh"
#include "Physica/Core/Utils/CUDA/device_obj.cuh"
#include "HostAllocator.h"

namespace Physica {
    template<class T>
    class DeviceAllocator {
        using This = DeviceAllocator<T>;
    public:
        using value_type = std::allocator_traits<This>::value_type;
        using pointer = std::allocator_traits<This>::pointer;
    public:
        DeviceAllocator() = default;
        DeviceAllocator(const This&) = default;
        DeviceAllocator(This&&) noexcept = default;
        ~DeviceAllocator() = default;
        /* Operators */
        This& operator=(const This&) noexcept = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        [[nodiscard, gnu::returns_nonnull]] __host__ __device__ pointer allocate(size_t n) noexcept;
        __host__ __device__ void deallocate(pointer p, size_t n) noexcept;
    };

    template<class T>
    __host__ __device__ auto DeviceAllocator<T>::allocate(size_t n) noexcept -> pointer {
        assert(n > 0 && "[Error]: Allocate nothing");
        size_t size = n * sizeof(value_type);
        void* p{};
        if constexpr (IsDevice())
            p = ::operator new(size);
        else {
            if constexpr (CUDADevAttr::MemoryPoolsSupported)
                cudaMallocAsync(&p, size, CUDAContext::getInstance());
            else
                cudaMalloc(&p, size);
        }
        return reinterpret_cast<pointer>(p);
    }

    template<class T>
    __host__ __device__ void DeviceAllocator<T>::deallocate(pointer p, [[maybe_unused]] size_t n) noexcept {
        if constexpr (IsDevice())
            ::operator delete(p);
        else {
            // No unnecessary cuda api call to make profiler output cleaner
            if (p != nullptr) {
                if constexpr (CUDADevAttr::MemoryPoolsSupported)
                    cudaFreeAsync(p, CUDAContext::getInstance());
                else
                    cudaFree(p);
            }
        }
    }
}

namespace std {
    template<class T>
    struct allocator_traits<Physica::DeviceAllocator<T>> {
        template<bool IsClass>
        struct ValueType { using Type = T::device_obj_type; };

        template<>
        struct ValueType<false> { using Type = T; };
    public:
        using allocator_type = Physica::DeviceAllocator<T>;
        using value_type = ValueType<std::is_class<T>::value>::Type;
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
    private:
        constexpr static bool IsTrivialC = std::is_trivially_constructible<value_type>::value;
        constexpr static bool IsTrivialD = std::is_trivially_destructible<value_type>::value;
        constexpr static bool IsTrivialCP = std::is_trivially_copyable<value_type>::value;
        static_assert((IsTrivialC && IsTrivialD) || IsTrivialCP, "[Error]: Must either be able to construct an object on device or copy it to device");
    public:
        [[nodiscard]] __host__ __device__ static pointer allocate(allocator_type& a, size_type n) {
            return a.allocate(n);
        }

        __host__ __device__ static void deallocate(allocator_type& a, pointer p, size_type n) noexcept {
            a.deallocate(p, n);
        }

        __host__ __device__ static void construct(allocator_type&, pointer p, auto&&... args) {
            using namespace Physica;
            if constexpr (IsDevice())
                ::new (static_cast<void*>(p)) value_type(std::forward<decltype(args)>(args)...);
            else if constexpr (sizeof...(args) != 0 || !IsTrivialC) {
                value_type temp(std::forward<decltype(args)>(args)...);
                check(cudaMemcpyAsync(p, &temp, sizeof(value_type), cudaMemcpyHostToDevice, CUDAContext::getInstance()));
                CUDAContext::getInstance().wait();
            }
        }

        __host__ __device__ static void destroy(allocator_type&, pointer p) {
            using namespace Physica;
            if constexpr (IsDevice())
                p->~value_type();
            else if constexpr (!IsTrivialD) {
                value_type temp{};
                cudaMemcpyAsync(&temp, p, sizeof(value_type), cudaMemcpyDeviceToHost, CUDAContext::getInstance());
                CUDAContext::getInstance().wait();
            }
        }

         __host__ __device__ static constexpr size_type max_size([[maybe_unused]] const allocator_type& a) noexcept {
            return std::numeric_limits<size_type>::max() / sizeof(value_type);
        }

         __host__ __device__ static allocator_type select_on_container_copy_construction(const allocator_type& a) {
            return a;
        }
    };
}
