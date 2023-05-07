/*
 * Copyright 2022-2023 WeiBo He.
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
#include "Physica/Utils/CUDA/DebugUtil.cuh"

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
    class device_obj<Array<T, Length, Capacity, Allocator>>
            : public Internal::ArrayBase<device_obj<Array<T, Length, Capacity, Allocator>>, DeviceAllocator<T>> {
        static_assert(Length != Dynamic, "[Error]: Dynamic length is not implemented");
        using host_obj = Array<T, Length, Capacity, Allocator>;
        using Base = Internal::ArrayBase<device_obj<Array<T, Length, Capacity, Allocator>>, DeviceAllocator<T>>;
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
    private:
        pointer d_data;
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
        device_obj& operator=(device_obj obj) noexcept;
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
        void toHost(host_obj& obj) const;
        template<class... Args> void resize([[maybe_unused]] size_t size, [[maybe_unused]] Args... args) { assert(size == Length); }
        void swap(device_obj& obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr static size_t size() { return Length; }
        [[nodiscard]] __host__ __device__ constexpr static size_t getLength() { return Length; }
        [[nodiscard]] __host__ __device__ constexpr static size_t getCapacity() { return Capacity; }
        using Base::empty;
        [[nodiscard]] __host__ __device__ pointer data() noexcept { return d_data; }
        [[nodiscard]] __host__ __device__ const_pointer data() const noexcept { return d_data; }
    };

    template<class T, size_t Length, size_t Capacity, class Allocator>
    device_obj<Array<T, Length, Capacity, Allocator>>::device_obj() : alloc() {
        d_data = alloc.allocate(Length);
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    template<class... Args>
    device_obj<Array<T, Length, Capacity, Allocator>>::device_obj(size_t length_, Args... args)
            : device_obj(host_obj(length_, std::forward<Args>(args)...).toDevice()) {}

    template<class T, size_t Length, size_t Capacity, class Allocator>
    device_obj<Array<T, Length, Capacity, Allocator>>::device_obj(const host_obj& array) : device_obj() {
        cudaError_t err;
        if constexpr (std::is_trivial<T>::value)
            err = cudaMemcpy(d_data, array.data(), Length * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyHostToDevice);
        else {
            Array<ValueType, Dynamic, Dynamic> buffer(Length);
            for (size_t i = 0; i < Length; ++i)
                buffer[i] = array[i].toDevice();
            err = cudaMemcpy(d_data, buffer.data(), Length * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyHostToDevice);
            buffer.get_allocator().deallocate(buffer.release(), Length);
        }
        if (err != cudaError_t::cudaSuccess)
            throw Core::CudaException(err);
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    device_obj<Array<T, Length, Capacity, Allocator>>::device_obj(const device_obj<Array<T, Length, Capacity, Allocator>>& obj)
            : alloc(obj.alloc) {
        d_data = alloc.allocate(Capacity);
        cudaError_t err;
        if constexpr (std::is_trivial<T>::value)
            err = cudaMemcpy(d_data, obj.d_data, Length * sizeof(T), cudaMemcpyKind::cudaMemcpyDeviceToDevice);
        else {
            Array<ValueType, Dynamic, Dynamic> buffer(Length);
            cudaMemcpy(buffer.data(), obj.d_data, Length * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyDeviceToHost);
            Array<ValueType, Dynamic, Dynamic> buffer1 = buffer;
            err = cudaMemcpy(d_data, buffer1.data(), Length * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyHostToDevice);
            buffer.get_allocator().deallocate(buffer.release(), Length);
            buffer1.get_allocator().deallocate(buffer1.release(), Length);
        }

        if (err != cudaError_t::cudaSuccess)
            throw Core::CudaException(err);
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    device_obj<Array<T, Length, Capacity, Allocator>>::device_obj(device_obj<Array<T, Length, Capacity, Allocator>>&& obj) noexcept
            : d_data(obj.d_data), alloc(std::move(obj.alloc)) {
        obj.d_data = nullptr;
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    device_obj<Array<T, Length, Capacity, Allocator>>::~device_obj() {
        if constexpr (!std::is_trivial<T>::value) {
            Array<ValueType, Dynamic, Dynamic> buffer(Length);
            cudaMemcpy(buffer.data(), d_data, Length * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyDeviceToHost);
        }
        alloc.deallocate(d_data, Length);
        d_data = nullptr;
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    device_obj<Array<T, Length, Capacity, Allocator>>&
    device_obj<Array<T, Length, Capacity, Allocator>>::operator=(device_obj<Array<T, Length, Capacity, Allocator>> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    typename device_obj<Array<T, Length, Capacity, Allocator>>::host_obj
    device_obj<Array<T, Length, Capacity, Allocator>>::toHost() const {
        host_obj result(getLength());
        toHost(result);
        return result;
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    void device_obj<Array<T, Length, Capacity, Allocator>>::toHost(host_obj& obj) const {
        cudaError_t err;
        if constexpr (std::is_trivial<T>::value)
            err = cudaMemcpy(obj.data(), d_data, getLength() * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyDeviceToHost);
        else {
            Array<ValueType, Dynamic, Dynamic> buffer(getLength());
            err = cudaMemcpy(buffer.data(), d_data, getLength() * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyDeviceToHost);
            for (size_t i = 0; i < getLength(); ++i)
                obj[i] = buffer[i].toHost();
        }
        if (err != cudaError_t::cudaSuccess)
            throw Core::CudaException(err);
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    void device_obj<Array<T, Length, Capacity, Allocator>>::swap(device_obj<Array<T, Length, Capacity, Allocator>>& obj) noexcept {
        std::swap(d_data, obj.d_data);
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    inline device_obj<Array<T, Length, Capacity, Allocator>> Array<T, Length, Capacity, Allocator>::toDevice() const {
        return device_obj<Array<T, Length, Capacity, Allocator>>(*this);
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    void Array<T, Length, Capacity, Allocator>::toDevice(device_obj<Array<T, Length, Capacity, Allocator>>& obj) const {
        cudaError_t err;
        if constexpr (std::is_trivial<T>::value)
            err = cudaMemcpy(obj.d_data, data(), Length * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyHostToDevice);
        else {
            Array<ValueType, Dynamic, Dynamic> buffer(Length);
            for (size_t i = 0; i < Length; ++i)
                buffer[i] = arr[i].toDevice();
            err = cudaMemcpy(obj.d_data, buffer.data(), Length * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyHostToDevice);
            buffer.get_allocator().deallocate(buffer.release(), Length);
        }
        if (err != cudaError_t::cudaSuccess)
            throw Core::CudaException(err);
    }

    template<class T, class Allocator>
    class device_obj<Array<T, Dynamic, Dynamic, Allocator>>
            : public Internal::ArrayBase<device_obj<Array<T, Dynamic, Dynamic, Allocator>>, DeviceAllocator<T>> {
        using host_obj = Array<T, Dynamic, Dynamic, Allocator>;
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
        device_obj& operator=(device_obj obj) noexcept;
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
        inline void toHost(host_obj& obj) const;
        void reserve(size_t size);
        template<class... Args> void resize(size_t size, Args... args);
        void swap(device_obj& obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t size() const noexcept { return getLength(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return length; }
        [[nodiscard]] __host__ __device__ size_t getCapacity() const noexcept { return capacity; }
        using Base::empty;
        [[nodiscard]] __host__ __device__ pointer data() noexcept { return d_data; }
        [[nodiscard]] __host__ __device__ const_pointer data() const noexcept { return d_data; }
    };

    template<class T, class Allocator>
    device_obj<Array<T, Dynamic, Dynamic, Allocator>>::device_obj() : d_data(nullptr), length(0), capacity(0) {}

    template<class T, class Allocator>
    template<class... Args>
    device_obj<Array<T, Dynamic, Dynamic, Allocator>>::device_obj(size_t length_, Args... args)
            : device_obj(host_obj(length_, std::forward<Args>(args)...).toDevice()) {}

    template<class T, class Allocator>
    device_obj<Array<T, Dynamic, Dynamic, Allocator>>::device_obj(const host_obj& array)
            : length(array.getLength()), capacity(array.getCapacity()) {
        d_data = alloc.allocate(capacity);
        cudaError_t err;
        if constexpr (std::is_trivial<T>::value)
            err = cudaMemcpy(d_data, array.data(), length * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyHostToDevice);
        else {
            Array<ValueType, Dynamic, Dynamic> buffer(length);
            for (size_t i = 0; i < length; ++i)
                buffer[i] = array[i].toDevice();
            err = cudaMemcpy(d_data, buffer.data(), length * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyHostToDevice);
            buffer.get_allocator().deallocate(buffer.release(), length);
        }
        if (err != cudaError_t::cudaSuccess)
            throw Core::CudaException(err);
    }

    template<class T, class Allocator>
    device_obj<Array<T, Dynamic, Dynamic, Allocator>>::device_obj(const device_obj<Array<T, Dynamic, Dynamic, Allocator>>& obj)
            : length(obj.getLength()), capacity(obj.getCapacity()), alloc(obj.alloc) {
        d_data = alloc.allocate(capacity);
        cudaError_t err;
        if constexpr (std::is_trivial<T>::value)
            err = cudaMemcpy(d_data, obj.d_data, length * sizeof(T), cudaMemcpyKind::cudaMemcpyDeviceToDevice);
        else {
            Array<ValueType, Dynamic, Dynamic> buffer(length);
            cudaMemcpy(buffer.data(), obj.d_data, length * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyDeviceToHost);
            Array<ValueType, Dynamic, Dynamic> buffer1 = buffer;
            err = cudaMemcpy(d_data, buffer1.data(), length * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyHostToDevice);
            buffer.get_allocator().deallocate(buffer.release(), length);
            buffer1.get_allocator().deallocate(buffer1.release(), length);
        }

        if (err != cudaError_t::cudaSuccess)
            throw Core::CudaException(err);
    }

    template<class T, class Allocator>
    device_obj<Array<T, Dynamic, Dynamic, Allocator>>::device_obj(device_obj<Array<T, Dynamic, Dynamic, Allocator>>&& obj) noexcept
            : d_data(obj.d_data), length(obj.length), capacity(obj.capacity), alloc(std::move(obj.alloc)) {
        obj.d_data = nullptr;
        obj.length = obj.capacity = 0;
    }

    template<class T, class Allocator>
    device_obj<Array<T, Dynamic, Dynamic, Allocator>>::~device_obj() {
        if constexpr (!std::is_trivial<T>::value) {
            Array<ValueType, Dynamic, Dynamic> buffer(length);
            cudaMemcpy(buffer.data(), d_data, length * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyDeviceToHost);
        }
        alloc.deallocate(d_data, length);
        d_data = nullptr;
        length = capacity = 0;
    }

    template<class T, class Allocator>
    device_obj<Array<T, Dynamic, Dynamic, Allocator>>&
    device_obj<Array<T, Dynamic, Dynamic, Allocator>>::operator=(device_obj<Array<T, Dynamic, Dynamic, Allocator>> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class T, class Allocator>
    Array<T, Dynamic, Dynamic, Allocator> device_obj<Array<T, Dynamic, Dynamic, Allocator>>::toHost() const {
        host_obj result(getLength());
        cudaError_t err;
        if constexpr (std::is_trivial<T>::value)
            err = cudaMemcpy(result.data(), d_data, getLength() * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyDeviceToHost);
        else {
            Array<ValueType, Dynamic, Dynamic> buffer(length);
            err = cudaMemcpy(buffer.data(), d_data, length * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyDeviceToHost);
            for (size_t i = 0; i < length; ++i)
                result[i] = buffer[i].toHost();
        }
        if (err != cudaError_t::cudaSuccess)
            throw Core::CudaException(err);
        return result;
    }

    template<class T, class Allocator>
    void device_obj<Array<T, Dynamic, Dynamic, Allocator>>::toHost(host_obj& obj) const {
        obj = toHost();
    }

    template<class T, class Allocator>
    void device_obj<Array<T, Dynamic, Dynamic, Allocator>>::reserve(size_t size) {
        assert(size > getCapacity());
        Array<ValueType, Dynamic, Dynamic> buffer(getLength());
        cudaCheck(cudaMemcpy(buffer.data(), d_data, getLength() * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyDeviceToHost));
        alloc.deallocate(d_data, capacity);
        d_data = alloc.allocate(size);
        cudaCheck(cudaMemcpy(d_data, buffer.data(), getLength() * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyHostToDevice));
        buffer.get_allocator().deallocate(buffer.release(), length);
        capacity = size;
    }

    template<class T, class Allocator>
    template<class... Args>
    void device_obj<Array<T, Dynamic, Dynamic, Allocator>>::resize(size_t size, Args... args) {
        if (capacity < size) {
            reserve(size);
        }

        Array<ValueType, Dynamic, Dynamic> buffer(std::max(length, size));
        cudaCheck(cudaMemcpy(buffer.data(), d_data, getLength() * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyDeviceToHost));
        if (length > size) {
            if constexpr (!std::is_trivial<T>::value)
                for (size_t i = size; i < length; ++i)
                    buffer[i].~ValueType();
            length = size;
        }
        else {
            for (; length < size; ++length)
                buffer.get_allocator().construct(buffer.data() + length, std::forward<Args>(args)...);
            cudaCheck(cudaMemcpy(d_data, buffer.data(), getLength() * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyHostToDevice));
        }
        buffer.get_allocator().deallocate(buffer.release(), length);
    }

    template<class T, class Allocator>
    void device_obj<Array<T, Dynamic, Dynamic, Allocator>>::swap(device_obj& obj) noexcept {
        std::swap(d_data, obj.d_data);
        std::swap(length, obj.length);
        std::swap(capacity, obj.capacity);
    }

    template<class T, class Allocator>
    inline device_obj<Array<T, Dynamic, Dynamic, Allocator>> Array<T, Dynamic, Dynamic, Allocator>::toDevice() const {
        return device_obj<Array<T, Dynamic, Dynamic, Allocator>>(*this);
    }

    template<class T, class Allocator>
    inline void Array<T, Dynamic, Dynamic, Allocator>::toDevice(device_obj<Array<T, Dynamic, Dynamic, Allocator>>& obj) const {
        obj = toDevice();
    }
}
