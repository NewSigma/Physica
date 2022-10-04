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
    namespace Internal {
        template<class T, size_t Length, size_t Capacity, class Allocator>
        class Traits<device_obj<Array<T, Length, Capacity, Allocator>>> {
        public:
            constexpr static size_t ArrayLength = Length;
            constexpr static size_t ArrayCapacity = Capacity;
            using AllocatorType = DeviceAllocator<T>;
            using ValueType = typename AllocatorType::value_type;
        };
    }

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

    template<class T, size_t Length, size_t Capacity, class Allocator_>
    inline device_obj<Array<T, Length, Capacity, Allocator_>> Array<T, Length, Capacity, Allocator_>::toDevice() const {
        return device_obj<Array<T, Length, Capacity, Allocator_>>(*this);
    }

    template<class T, class Allocator_>
    class device_obj<Array<T, Dynamic, Dynamic, Allocator_>>
            : private Internal::ArrayBase<device_obj<Array<T, Dynamic, Dynamic, Allocator_>>, DeviceAllocator<T>> {
        using host_obj = Array<T, Dynamic, Dynamic, Allocator_>;
        using This = device_obj<host_obj>;
        using Base = Internal::ArrayBase<device_obj<host_obj>, DeviceAllocator<T>>;
        using typename Base::allocator_type;
        using typename Base::ValueType;
        using typename Base::pointer;
        using typename Base::const_pointer;
        using typename Base::lvalue_reference;
        using typename Base::const_lvalue_reference;
        using typename Base::rvalue_reference;
        using typename Base::Iterator;
        using typename Base::ConstIterator;
        using typename Base::ReverseIterator;
        using typename Base::ConstReverseIterator;

        pointer d_data;
        size_t length;
        size_t capacity;
        allocator_type alloc;
    public:
        device_obj();
        template<class... Args>
        explicit device_obj(size_t length_, Args... args);
        explicit device_obj(const host_obj& array);
        device_obj(const device_obj& obj);
        device_obj(device_obj&& obj) noexcept;
        ~device_obj();
        /* Operators */
        device_obj& operator=(device_obj other) noexcept;
        [[nodiscard]] __device__ lvalue_reference operator[](size_t index) { return Base::operator[](index); }
        [[nodiscard]] __device__ const_lvalue_reference operator[](size_t index) const { return Base::operator[](index); }
        /* Iterator */
        __device__ Iterator begin() noexcept { return Base::begin(); }
        __device__ ConstIterator begin() const noexcept { return Base::cbegin(); }
        __device__ ConstIterator cbegin() const noexcept { return Base::cbegin(); }
        __device__ Iterator end() noexcept { return Base::end(); }
        __device__ ConstIterator end() const noexcept { return Base::cend(); }
        __device__ ConstIterator cend() const noexcept { return Base::cend(); }
        __device__ ReverseIterator rbegin() noexcept { return Base::rbegin(); }
        __device__ ConstReverseIterator rbegin() const noexcept { return Base::crbegin(); }
        __device__ ConstReverseIterator crbegin() const noexcept { return Base::cbegin(); }
        __device__ ReverseIterator rend() noexcept { return Base::rend(); }
        __device__ ConstReverseIterator rend() const noexcept { return Base::crend(); }
        __device__ ConstReverseIterator crend() const noexcept { return Base::crend(); }
        /* Operations */
        [[nodiscard]] host_obj toHost() const;
        void swap(device_obj& obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t size() const noexcept { return getLength(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return length; }
        [[nodiscard]] __host__ __device__ size_t getCapacity() const noexcept { return capacity; }
        using Base::empty;
        [[nodiscard]] __host__ __device__ pointer data() noexcept { return d_data; }
    };

    template<class T, class Allocator_>
    device_obj<Array<T, Dynamic, Dynamic, Allocator_>>::device_obj() : d_data(nullptr), length(0), capacity(0) {}

    template<class T, class Allocator_>
    template<class... Args>
    device_obj<Array<T, Dynamic, Dynamic, Allocator_>>::device_obj(size_t length_, Args... args)
            : device_obj(host_obj(length_, std::forward<Args>(args)...).toDevice()) {}

    template<class T, class Allocator_>
    device_obj<Array<T, Dynamic, Dynamic, Allocator_>>::device_obj(const host_obj& array)
            : length(array.getLength()), capacity(array.getCapacity()) {
        d_data = alloc.allocate(capacity);
        cudaError_t err;
        if constexpr (std::is_trivial<T>::value)
            err = cudaMemcpy(d_data, array.data(), length * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyHostToDevice);
        else {
            Array<ValueType, Dynamic, Dynamic, HostAllocator<ValueType>> buffer(length);
            for (size_t i = 0; i < length; ++i)
                buffer[i] = array[i].toDevice();
            err = cudaMemcpy(d_data, buffer.data(), length * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyHostToDevice);
            buffer.get_allocator().deallocate(buffer.release(), length);
        }
        if (err != cudaError_t::cudaSuccess)
            throw Core::CudaException(err);
    }

    template<class T, class Allocator_>
    device_obj<Array<T, Dynamic, Dynamic, Allocator_>>::device_obj(const device_obj<Array<T, Dynamic, Dynamic, Allocator_>>& obj)
            : length(obj.getLength()), capacity(obj.getCapacity()), alloc(obj.alloc) {
        d_data = alloc.allocate(capacity);
        const auto code = cudaMemcpy(d_data, obj.d_data, length * sizeof(T), cudaMemcpyKind::cudaMemcpyDeviceToDevice);
        if (code != cudaError_t::cudaSuccess)
            throw Core::CudaException(code);
    }

    template<class T, class Allocator_>
    device_obj<Array<T, Dynamic, Dynamic, Allocator_>>::device_obj(device_obj<Array<T, Dynamic, Dynamic, Allocator_>>&& obj) noexcept
            : d_data(obj.d_data), length(obj.length), capacity(obj.capacity), alloc(std::move(obj.alloc)) {
        obj.d_data = nullptr;
        obj.length = obj.capacity = 0;
    }

    template<class T, class Allocator_>
    device_obj<Array<T, Dynamic, Dynamic, Allocator_>>::~device_obj() {
        if constexpr (!std::is_trivial<T>::value) {
            Array<ValueType, Dynamic, Dynamic, HostAllocator<ValueType>> buffer(length);
            cudaMemcpy(buffer.data(), d_data, length * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyDeviceToHost);
        }
        alloc.deallocate(d_data, length);
        d_data = nullptr;
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
        cudaError_t err;
        if constexpr (std::is_trivial<T>::value)
            err = cudaMemcpy(result.data(), d_data, getLength() * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyDeviceToHost);
        else {
            Array<ValueType, Dynamic, Dynamic, HostAllocator<ValueType>> buffer(length);
            err = cudaMemcpy(buffer.data(), d_data, length * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyDeviceToHost);
            for (size_t i = 0; i < length; ++i)
                result[i] = buffer[i].toHost();
        }
        if (err != cudaError_t::cudaSuccess)
            throw Core::CudaException(err);
        return result;
    }

    template<class T, class Allocator_>
    void device_obj<Array<T, Dynamic, Dynamic, Allocator_>>::swap(device_obj& obj) noexcept {
        std::swap(d_data, obj.d_data);
        std::swap(length, obj.length);
        std::swap(capacity, obj.capacity);
    }

    template<class T, class Allocator_>
    inline device_obj<Array<T, Dynamic, Dynamic, Allocator_>> Array<T, Dynamic, Dynamic, Allocator_>::toDevice() const {
        return device_obj<Array<T, Dynamic, Dynamic, Allocator_>>(*this);
    }
}
