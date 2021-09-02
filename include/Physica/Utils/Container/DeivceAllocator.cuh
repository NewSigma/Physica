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
#include <memory>
#include <thrust/device_ptr.h>
#include <cuda_runtime.h>

namespace Physica::Utils {
    template<class T> class DeviceAllocator;
}

namespace std {
    template<class T>
    struct allocator_traits<Physica::Utils::DeviceAllocator<T>> {
    public:
        using allocator_type = Physica::Utils::DeviceAllocator<T>;
        using value_type = T;
        using pointer = thrust::device_ptr<T>;
        using const_pointer = const thrust::device_ptr<T>;
        using void_pointer = thrust::device_ptr<void>;
        using const_void_pointer = const thrust::device_ptr<void>;
        using lvalue_reference = thrust::device_reference<T>;
        using const_lvalue_reference = const T&;
        using rvalue_reference = T&&;
        using size_type = size_t;
        using difference_type = std::ptrdiff_t;
        using propagate_on_container_copy_assignment = std::false_type;
        using propagate_on_container_move_assignment = std::false_type;
        using propagate_on_container_swap = std::false_type;
        using is_always_equal = typename std::is_empty<allocator_type>::type;
        template<class U>
        using rebind_alloc = Physica::Utils::DeviceAllocator<U>;
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
    class DeviceAllocator {
    public:
        using value_type = T;
        using pointer = typename std::allocator_traits<DeviceAllocator>::pointer;
        using size_type = typename std::allocator_traits<DeviceAllocator>::size_type;
        using difference_type = typename std::allocator_traits<DeviceAllocator>::difference_type;
        using propagate_on_container_move_assignment = std::true_type;
    public:
        DeviceAllocator() noexcept = default;
        DeviceAllocator(const DeviceAllocator&) noexcept = default;
        DeviceAllocator(DeviceAllocator&&) noexcept = delete;
        ~DeviceAllocator() = default;
        /* Operators */
        DeviceAllocator& operator=(const DeviceAllocator&) noexcept = default;
        DeviceAllocator& operator=(DeviceAllocator&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __host__ __device__ static pointer allocate(size_t n);
        __host__ __device__ static void deallocate(pointer p, size_t n) noexcept;
        template<class... Args>
        __host__ __device__ void construct(pointer p, Args&&... args);
        __host__ __device__ void destroy(pointer p);
    };

    template<class T>
    __host__ __device__ typename DeviceAllocator<T>::pointer DeviceAllocator<T>::allocate(size_t n) {
    #ifdef __CUDA__ARCH__
        auto* p = reinterpret_cast<T*>(malloc(n * sizeof(T)));
        return thrust::device_ptr(p);
    #else
        T* p;
        cudaMalloc(&p, n * sizeof(T));
        return pointer(p);
    #endif
    }

    template<class T>
    __host__ __device__ void DeviceAllocator<T>::deallocate(pointer p, [[maybe_unused]] size_t n) noexcept {
    #ifdef __CUDA__ARCH__
        free(p.get());
    #else
        cudaFree(p.get());
    #endif
    }

    template<class T>
    template<class... Args>
    __host__ __device__ void DeviceAllocator<T>::construct(pointer p, Args&&... args) {
    #ifdef __CUDA__ARCH__
        ::new (static_cast<void*>(p.get())) T(std::forward<Args>(args)...);
    #else
        T temp(args...);
        cudaMemcpy(p.get(), &temp, sizeof(T), cudaMemcpyHostToDevice);
        if constexpr (!std::is_trivial<T>::value)
            temp.release(); //Ownership has been given to device
    #endif
    }

    template<class T>
    __host__ __device__ void DeviceAllocator<T>::destroy(pointer p) {
    #ifdef __CUDA__ARCH__
        p->~T();
    #else
        T temp;
        cudaMemcpy(&temp, p.get(), sizeof(T), cudaMemcpyDeviceToHost);
    #endif
    }
}
