/*
 * Copyright 2022-2024 Weibo He.
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

#include <Physica/Utils/Allocator/DeviceAllocator.cuh>
#include <Physica/Utils/CUDA/PlainStruct.h>
#include <Physica/Core/Exception/CudaException.cuh>

namespace Physica::Utils {
    template<class T, size_t Length, size_t Capacity, class Allocator>
    class device_obj<Array<T, Length, Capacity, Allocator>> : public Array<T, Length, Capacity, Allocator> {
        static_assert(Length != Dynamic, "[Error]: Dynamic length is not implemented");
        static_assert(std::is_trivial<T>::value, "[Error]: Fixed size array with non-trivial element is not supported, it is seldom used on cuda");
        using host_obj = Array<T, Length, Capacity, Allocator>;
        using Base = host_obj;
    public:
        using host_obj::host_obj;
        device_obj() = default;
        __host__ __device__ device_obj(const host_obj& obj) : host_obj(obj) {}
        device_obj(const device_obj& obj) = default;
        device_obj(device_obj&& obj) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        __host__ __device__ device_obj& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        /* Operations */
        using Base::resize;
        [[nodiscard]] __host__ __device__ host_obj toHost() const { return *this; }
        __host__ __device__ void toHost(host_obj& obj) const { obj = *this; }
        __host__ __device__ void swap(device_obj& __restrict obj) noexcept { host_obj::swap(obj); }
    };

    template<class T, size_t Length, size_t Capacity, class Allocator>
    inline device_obj<Array<T, Length, Capacity, Allocator>> Array<T, Length, Capacity, Allocator>::toDevice() const {
        return device_obj<Array<T, Length, Capacity, Allocator>>(*this);
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    inline void Array<T, Length, Capacity, Allocator>::toDevice(device_obj<Array<T, Length, Capacity, Allocator>>& obj) const {
        obj = *this;
    }

    template<class T, class Allocator>
    class device_obj<Array<T, Dynamic, Dynamic, Allocator>>
            : public ArrayBase<device_obj<Array<T, Dynamic, Dynamic, Allocator>>, DeviceAllocator<T>> {
        static_assert(!is_device_obj<T>::value, "[Error]: Nested device_obj is not allowed");
        using host_obj = Array<T, Dynamic, Dynamic, Allocator>;
        using This = device_obj<host_obj>;
        using Base = ArrayBase<device_obj<host_obj>, DeviceAllocator<T>>;
    public:
        constexpr static bool isTrivial = std::is_trivial<T>::value;
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
        using PlainElemType = typename std::conditional<isTrivial, ValueType, Physica::PlainStruct<ValueType>>::type;
        using PlainElemAllocator = typename ChangeAllocatorValueType<Allocator, PlainElemType>::Type;
        using PlainHostObj = Array<PlainElemType, Dynamic, Dynamic, PlainElemAllocator>;
    private:
        pointer d_data;
        size_t length;
        size_t capacity;
        allocator_type alloc;
    public:
        device_obj();
        template<class... Args>
        explicit device_obj(size_t length_, Args&&... args);
        explicit device_obj(const host_obj& array);
        device_obj(const device_obj& obj);
        device_obj(device_obj&& obj) noexcept;
        ~device_obj();
        /* Operators */
        device_obj& operator=(device_obj obj) noexcept;
        [[nodiscard]] __device__ lvalue_reference operator[](size_t index) { return Base::operator[](index); }
        [[nodiscard]] __device__ const_lvalue_reference operator[](size_t index) const { return Base::operator[](index); }
        /* Iterators */
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
        [[nodiscard]] PlainHostObj toPlainHost() const;
        [[nodiscard]] host_obj toHost() const;
        inline void toHost(host_obj& obj) const;
        void reserve(size_t size);
        template<class... Args> void resize(size_t size, Args&&... args);
        void swap(device_obj& __restrict obj) noexcept;
        [[nodiscard]] pointer release() noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t size() const noexcept { return getLength(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return length; }
        [[nodiscard]] __host__ __device__ size_t getCapacity() const noexcept { return capacity; }
        using Base::empty;
        [[nodiscard]] __host__ __device__ pointer data() noexcept { return d_data; }
        [[nodiscard]] __host__ __device__ const_pointer data() const noexcept { return d_data; }
    private:
        using Base::operator[];
        /* Iterators */
        using Base::begin;
        using Base::cbegin;
        using Base::end;
        using Base::cend;
        using Base::rbegin;
        using Base::crbegin;
        using Base::rend;
        using Base::crend;
    };

    template<class T, class Allocator>
    device_obj<Array<T, Dynamic, Dynamic, Allocator>>::device_obj() : d_data(nullptr), length(0), capacity(0) {}

    template<class T, class Allocator>
    template<class... Args>
    device_obj<Array<T, Dynamic, Dynamic, Allocator>>::device_obj(size_t length_, Args&&... args)
            : device_obj(host_obj(length_, std::forward<Args>(args)...).toDevice()) {}

    template<class T, class Allocator>
    device_obj<Array<T, Dynamic, Dynamic, Allocator>>::device_obj(const host_obj& array)
            : length(array.getLength()), capacity(array.getCapacity()) {
        d_data = alloc.allocate(capacity);
        array.toDevice(*this);
    }
    /**
     * Do not launch kernel before copy is finished, so memcpy is async.
     */
    template<class T, class Allocator>
    device_obj<Array<T, Dynamic, Dynamic, Allocator>>::device_obj(const device_obj<Array<T, Dynamic, Dynamic, Allocator>>& obj)
            : length(obj.getLength()), capacity(obj.getCapacity()), alloc(obj.alloc) {
        d_data = alloc.allocate(capacity);
        auto& stream = Core::StreamPool::getStream();
        if constexpr (isTrivial)
            cudaMemcpyAsync(d_data, obj.d_data, length * sizeof(T), cudaMemcpyKind::cudaMemcpyDeviceToDevice, stream);
        else {
            Array<ValueType, Dynamic, Dynamic> buffer(length);
            cudaMemcpyAsync(buffer.data(), obj.d_data, length * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyDeviceToHost, stream);
            stream.wait();
            Array<ValueType, Dynamic, Dynamic> buffer1 = buffer;
            cudaMemcpyAsync(d_data, buffer1.data(), length * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyHostToDevice, stream);
            buffer.get_allocator().deallocate(buffer.release(), length);
            stream.wait();
            buffer1.get_allocator().deallocate(buffer1.release(), length);
        }
    }

    template<class T, class Allocator>
    device_obj<Array<T, Dynamic, Dynamic, Allocator>>::device_obj(device_obj<Array<T, Dynamic, Dynamic, Allocator>>&& obj) noexcept
            : d_data(obj.d_data), length(obj.length), capacity(obj.capacity), alloc(std::move(obj.alloc)) {
        obj.d_data = nullptr;
        obj.length = obj.capacity = 0;
    }

    template<class T, class Allocator>
    device_obj<Array<T, Dynamic, Dynamic, Allocator>>::~device_obj() {
        if constexpr (!isTrivial) {
            Array<ValueType, Dynamic, Dynamic> buffer(length);
            auto& stream = Core::StreamPool::getStream();
            cudaMemcpyAsync(buffer.data(), d_data, length * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyDeviceToHost, stream);
            stream.wait();
        }
        alloc.deallocate(d_data, capacity);
        d_data = nullptr;
        length = capacity = 0;
    }

    template<class T, class Allocator>
    device_obj<Array<T, Dynamic, Dynamic, Allocator>>&
    device_obj<Array<T, Dynamic, Dynamic, Allocator>>::operator=(device_obj<Array<T, Dynamic, Dynamic, Allocator>> obj) noexcept {
        swap(obj);
        return *this;
    }
    /**
     * Synchronization is expected before this, copy is async to task stream
     */
    template<class T, class Allocator>
    typename device_obj<Array<T, Dynamic, Dynamic, Allocator>>::PlainHostObj
    device_obj<Array<T, Dynamic, Dynamic, Allocator>>::toPlainHost() const {
        PlainHostObj result(getLength());
        cudaCheck(cudaMemcpy(result.data(), (void*)d_data, getLength() * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyDeviceToHost));
        return result;
    }
    /**
     * Synchronization is expected before this, copy is async to task stream
     */
    template<class T, class Allocator>
    Array<T, Dynamic, Dynamic, Allocator> device_obj<Array<T, Dynamic, Dynamic, Allocator>>::toHost() const {
        host_obj result(length);
        if constexpr (isTrivial)
            cudaCheck(cudaMemcpy(result.data(), d_data, length * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyDeviceToHost));
        else {
            const auto buffer = toPlainHost();
            for (size_t i = 0; i < length; ++i)
                result[i] = buffer[i].getDerived().toHost();
        }
        return result;
    }

    template<class T, class Allocator>
    inline void device_obj<Array<T, Dynamic, Dynamic, Allocator>>::toHost(host_obj& obj) const {
        obj = toHost();
    }
    /**
     * Synchronization is expected before this, copy is async to task stream
     */
    template<class T, class Allocator>
    void device_obj<Array<T, Dynamic, Dynamic, Allocator>>::reserve(size_t size) {
        assert(size > getCapacity());
        const auto buffer = toPlainHost();
        alloc.deallocate(d_data, capacity);
        d_data = alloc.allocate(size);
        cudaCheck(cudaMemcpy(d_data, buffer.data(), getLength() * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyHostToDevice));
        capacity = size;
        cudaCheck(cudaStreamSynchronize(nullptr));
    }
    /**
     * Synchronization is expected before this, copy is async to task stream
     */
    template<class T, class Allocator>
    template<class... Args>
    void device_obj<Array<T, Dynamic, Dynamic, Allocator>>::resize(size_t size, Args&&... args) {
        if (size == length)
            return;
        if (capacity < size)
            reserve(size);

        if (length > size) {
            if constexpr (!isTrivial) {
                const size_t delta = length - size;
                Array<ValueType, Dynamic, Dynamic> buffer{};
                buffer.reserve(delta);
                cudaCheck(cudaMemcpy(buffer.data(), d_data + size, delta * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyDeviceToHost));
                buffer.setLength(delta);
            }
        }
        else {
            const size_t delta = size - length;
            Array<ValueType, Dynamic, Dynamic> buffer(delta);
            for (size_t i = 0; i < delta; ++i)
                buffer.get_allocator().construct(buffer.data() + i, std::forward<Args>(args)...);
            cudaCheck(cudaMemcpy(d_data + length, buffer.data(), delta * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyHostToDevice));
            buffer.get_allocator().deallocate(buffer.release(), delta);
            cudaCheck(cudaStreamSynchronize(nullptr));
        }
        length = size;
    }

    template<class T, class Allocator>
    void device_obj<Array<T, Dynamic, Dynamic, Allocator>>::swap(device_obj& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(d_data, obj.d_data);
        std::swap(length, obj.length);
        std::swap(capacity, obj.capacity);
    }

    template<class T, class Allocator>
    typename device_obj<Array<T, Dynamic, Dynamic, Allocator>>::pointer
    device_obj<Array<T, Dynamic, Dynamic, Allocator>>::release() noexcept {
        auto* copy = d_data;
        d_data = nullptr;
        length = capacity = 0;
        return copy;
    }

    template<class T, class Allocator>
    inline device_obj<Array<T, Dynamic, Dynamic, Allocator>> Array<T, Dynamic, Dynamic, Allocator>::toDevice() const {
        return device_obj<Array<T, Dynamic, Dynamic, Allocator>>(*this);
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Dynamic, Allocator>::toDevice(device_obj<Array<T, Dynamic, Dynamic, Allocator>>& obj) const {
        using ValueType = typename device_obj<This>::ValueType;
        constexpr bool isTrivial = device_obj<This>::isTrivial;
        const size_t length = getLength();
        obj.resize(length);

        auto& stream = Core::StreamPool::getStream();
        if constexpr (isTrivial) {
            cudaMemcpyAsync(obj.data(), this->data(), length * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyHostToDevice, stream);
            stream.wait();
        }
        else {
            Array<ValueType, Dynamic, Dynamic> buffer(length);
            for (size_t i = 0; i < length; ++i)
                buffer[i] = this->operator[](i).toDevice();
            cudaMemcpyAsync(obj.data(), buffer.data(), length * sizeof(ValueType), cudaMemcpyKind::cudaMemcpyHostToDevice, stream);
            stream.wait();
            buffer.get_allocator().deallocate(buffer.release(), length);
        }
    }
}

namespace Physica {
    using namespace Utils;

    template<class T, size_t Length, size_t Capacity, class Allocator>
    class Traits<Utils::device_obj<Array<T, Length, Capacity, Allocator>>> {
    public:
        constexpr static size_t ArrayLength = Length;
        constexpr static size_t ArrayCapacity = Capacity;
        using AllocatorType = DeviceAllocator<T>;
        using ValueType = typename AllocatorType::value_type;
    };
}
