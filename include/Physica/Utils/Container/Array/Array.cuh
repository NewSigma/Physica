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
    class device_obj<Array<T, Dynamic, Dynamic, Allocator_>> : public Internal::ArrayBase<device_obj<Array<T, Dynamic, Dynamic, Allocator_>>, DeviceAllocator<T>> {
        using host_obj = Array<T, Dynamic, Dynamic, Allocator_>;
        using This = device_obj<host_obj>;
        using Base = Internal::ArrayBase<device_obj<host_obj>, DeviceAllocator<T>>;
        using typename Base::allocator_type;

        thrust::device_ptr<T> data;
        size_t length;
        size_t capacity;
        allocator_type alloc;
    public:
        explicit device_obj(size_t length_);
        explicit device_obj(const host_obj& array);
        device_obj(const device_obj& obj);
        device_obj(device_obj&& obj) noexcept;
        ~device_obj();
        /* Operators */
        device_obj& operator=(device_obj other) noexcept;
        /* Operations */
        [[nodiscard]] host_obj toHost() const;
        void swap(device_obj& obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t size() const noexcept { return getLength(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return length; }
        [[nodiscard]] __host__ __device__ size_t getCapacity() const noexcept { return capacity; }
    };

    template<class T, class Allocator_>
    device_obj<Array<T, Dynamic, Dynamic, Allocator_>>::device_obj(size_t length_) : length(length_), capacity(length_) {
        data = alloc.allocate(capacity);
    }

    template<class T, class Allocator_>
    device_obj<Array<T, Dynamic, Dynamic, Allocator_>>::device_obj(const host_obj& array)
            : length(array.getLength()), capacity(array.getCapacity()) {
        data = alloc.allocate(capacity);
        const auto code = cudaMemcpy(data.get(), array.data(), length * sizeof(T), cudaMemcpyKind::cudaMemcpyHostToDevice);
        if (code != cudaError_t::cudaSuccess)
            throw Core::CudaException(code);
    }

    template<class T, class Allocator_>
    device_obj<Array<T, Dynamic, Dynamic, Allocator_>>::device_obj(const device_obj<Array<T, Dynamic, Dynamic, Allocator_>>& obj)
            : length(obj.getLength()), capacity(obj.getCapacity()), alloc(obj.alloc) {
        data = alloc.allocate(capacity);
        const auto code = cudaMemcpy(data.get(), obj.data.get(), length * sizeof(T), cudaMemcpyKind::cudaMemcpyDeviceToDevice);
        if (code != cudaError_t::cudaSuccess)
            throw Core::CudaException(code);
    }

    template<class T, class Allocator_>
    device_obj<Array<T, Dynamic, Dynamic, Allocator_>>::device_obj(device_obj<Array<T, Dynamic, Dynamic, Allocator_>>&& obj) noexcept
            : data(obj.data), length(obj.length), capacity(obj.capacity), alloc(std::move(obj.alloc)) {
        obj.data = nullptr;
        obj.length = obj.capacity = 0;
    }

    template<class T, class Allocator_>
    device_obj<Array<T, Dynamic, Dynamic, Allocator_>>::~device_obj() {
        alloc.deallocate(data, length);
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
        const auto code = cudaMemcpy(result.data(), data.get(), getLength() * sizeof(T), cudaMemcpyKind::cudaMemcpyDeviceToHost);
        if (code != cudaError_t::cudaSuccess)
            throw Core::CudaException(code);
        return result;
    }

    template<class T, class Allocator_>
    void device_obj<Array<T, Dynamic, Dynamic, Allocator_>>::swap(device_obj& obj) noexcept {
        std::swap(data, obj.data);
        std::swap(length, obj.length);
        std::swap(capacity, obj.capacity);
    }

    template<class T, size_t Length, size_t Capacity, class Allocator_>
    inline device_obj<Array<T, Length, Capacity, Allocator_>> Array<T, Length, Capacity, Allocator_>::toDevice() const {
        return device_obj<Array<T, Length, Capacity, Allocator_>>(*this);
    }

    template<class T, class Allocator_>
    inline device_obj<Array<T, Dynamic, Dynamic, Allocator_>> Array<T, Dynamic, Dynamic, Allocator_>::toDevice() const {
        return device_obj<Array<T, Dynamic, Dynamic, Allocator_>>(*this);
    }
}
