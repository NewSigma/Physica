/*
 * Copyright 2022 WeiBo He.
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

#include <thrust/device_ptr.h>
#include <cuda_runtime.h>
#include "Physica/Utils/Container/DeviceAllocator.cuh"

namespace Physica::Utils {
    template<class T, size_t Length, size_t Capacity, class Allocator>
    class device_obj<Array<T, Length, Capacity, Allocator>> : private Array<T, Length, Capacity, Allocator> {
        using host_obj = Array<T, Length, Capacity, Allocator>;
    public:
        /* Operators */
        using host_obj::operator=;
        /* Operations */
        using host_obj::swap;
        host_obj toHost() const { return host_obj(*this); }
        /* Getters */
        using host_obj::size;
        using host_obj::getLength;
        using host_obj::getCapacity;
    };

    template<class T, class Allocator_>
    class device_obj<Array<T, Dynamic, Dynamic, Allocator_>> {
        using This = device_obj<Array<T, Dynamic, Dynamic, Allocator_>>;
        using host_obj = Array<T, Dynamic, Dynamic, Allocator_>;
        using Allocator = DeviceAllocator<T>;

        thrust::device_ptr<T> data;
        size_t length;
        size_t capacity;
        Allocator alloc;
    public:
        device_obj(const host_obj& array);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj();
        /* Operators */
        device_obj& operator=(device_obj other) noexcept;
        host_obj toHost() const;
        void swap(device_obj& obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t size() { return getLength(); }
        [[nodiscard]] __host__ __device__ size_t getLength() { return length; }
        [[nodiscard]] __host__ __device__ size_t getCapacity() { return capacity; }
    };

    template<class T, class Allocator_>
    device_obj<Array<T, Dynamic, Dynamic, Allocator_>>::device_obj(const host_obj& array)
            : length(array.getLength()), capacity(array.getCapacity()) {
        data = Allocator::allocate(alloc, capacity);
        cudaMemcpy(data.get(), array.data(), length, cudaMemcpyKind::cudaMemcpyHostToDevice);
    }

    template<class T, class Allocator_>
    device_obj<Array<T, Dynamic, Dynamic, Allocator_>>::~device_obj() {
        Allocator::deallocate(data.get(), length);
        data = nullptr;
        length = capacity = 0;
    }

    template<class T, class Allocator_>
    device_obj<Array<T, Dynamic, Dynamic, Allocator_>>&
    device_obj<Array<T, Dynamic, Dynamic, Allocator_>>::operator=(device_obj<Array<T, Dynamic, Dynamic, Allocator_>> other) noexcept {
        swap(other);
        return *this;
    }

    template<class T, class Allocator_>
    Array<T, Dynamic, Dynamic, Allocator_> device_obj<Array<T, Dynamic, Dynamic, Allocator_>>::toHost() const {
        host_obj result(getCapacity());
        cudaMemcpy(result.data(), data.get(), getLength(), cudaMemcpyKind::cudaMemcpyDeviceToHost);
    }

    template<class T, class Allocator_>
    void device_obj<Array<T, Dynamic, Dynamic, Allocator_>>::swap(device_obj& obj) noexcept {
        std::swap(data, obj.data);
        std::swap(length, obj.length);
        std::swap(capacity, obj.capacity);
    }

    template<class T, size_t Length, size_t Capacity, class Allocator_>
    inline device_obj<Array<T, Length, Capacity, Allocator_>> Array<T, Length, Capacity, Allocator_>::toDevice() {
        return device_obj<Array<T, Length, Capacity, Allocator_>>(*this);
    }

    template<class T, class Allocator_>
    inline device_obj<Array<T, Dynamic, Dynamic, Allocator_>> Array<T, Dynamic, Dynamic, Allocator_>::toDevice() {
        return device_obj<Array<T, Dynamic, Dynamic, Allocator_>>(*this);
    }
}
