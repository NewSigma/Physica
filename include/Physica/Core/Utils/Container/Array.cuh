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

#include <Physica/PlainStruct.h>
#include <Physica/Core/Utils/Allocator/DeviceAllocator.cuh>
#include <Physica/Core/Exception/CUDA/CUDA.cuh>

namespace Physica::Core {
    template<class T, size_t Length, class Allocator>
    class device_obj<Array<T, Length, Allocator>> : public Array<T, Length, Allocator> {
        static_assert(Length != Dynamic, "[Error]: Dynamic length is not implemented");
        static_assert(std::is_trivial<T>::value, "[Error]: Fixed size array with non-trivial element is not supported, it is seldom used on cuda");
        using host_obj = Array<T, Length, Allocator>;
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

    template<class T, size_t Length, class Allocator>
    inline auto Array<T, Length, Allocator>::toDevice() const {
        return device_obj<This>(*this);
    }

    template<class T, size_t Length, class Allocator>
    inline auto Array<T, Length, Allocator>::toDeviceAsync() const {
        return toDevice();
    }

    template<class T, size_t Length, class Allocator>
    inline void Array<T, Length, Allocator>::toDevice(device_obj<This>& obj) const {
        obj = *this;
    }

    template<class T, size_t Length, class Allocator>
    inline void Array<T, Length, Allocator>::toDeviceAsync(device_obj<This>& obj) const {
        toDevice(obj);
    }

    template<class T, class Allocator>
    class device_obj<Array<T, Dynamic, Allocator>> : public ArrayBase<device_obj<Array<T, Dynamic, Allocator>>, DeviceAllocator<T>> {
        static_assert(!is_device_obj<T>::value, "[Error]: Nested device_obj is not allowed");
        using host_obj = Array<T, Dynamic, Allocator>;
        using This = device_obj<host_obj>;
        using Base = ArrayBase<device_obj<host_obj>, DeviceAllocator<T>>;
    public:
        constexpr static bool isTrivial = std::is_trivial<T>::value;
        using typename Base::allocator_type;
        using typename Base::value_type;
        using typename Base::pointer;
        using typename Base::const_pointer;
        using typename Base::lvalue_reference;
        using typename Base::const_lvalue_reference;
        using typename Base::rvalue_reference;

        using typename Base::ElemType;
        using PlainElemType = typename std::conditional<isTrivial, ElemType, Physica::PlainStruct<ElemType>>::type;
        using PlainElemAllocator = typename Allocator::template rebind_alloc<PlainElemType>;
        using PlainHostObj = Array<PlainElemType, Dynamic, PlainElemAllocator>;
    protected:
        using typename Base::FIteType;
        using typename Base::RIteType;
        using typename Base::CFIteType;
        using typename Base::CRIteType;
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
        __device__ FIteType begin() noexcept { return Base::begin(); }
        __device__ CFIteType begin() const noexcept { return Base::cbegin(); }
        __device__ CFIteType cbegin() const noexcept { return Base::cbegin(); }
        __device__ FIteType end() noexcept { return Base::end(); }
        __device__ CFIteType end() const noexcept { return Base::cend(); }
        __device__ CFIteType cend() const noexcept { return Base::cend(); }
        __device__ RIteType rbegin() noexcept { return Base::rbegin(); }
        __device__ CRIteType rbegin() const noexcept { return Base::crbegin(); }
        __device__ CRIteType crbegin() const noexcept { return Base::cbegin(); }
        __device__ RIteType rend() noexcept { return Base::rend(); }
        __device__ CRIteType rend() const noexcept { return Base::crend(); }
        __device__ CRIteType crend() const noexcept { return Base::crend(); }
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
    device_obj<Array<T, Dynamic, Allocator>>::device_obj() : d_data(nullptr), length(0), capacity(0) {}

    template<class T, class Allocator>
    template<class... Args>
    device_obj<Array<T, Dynamic, Allocator>>::device_obj(size_t length_, Args&&... args)
            : device_obj(host_obj(length_, std::forward<Args>(args)...).toDevice()) {}

    template<class T, class Allocator>
    device_obj<Array<T, Dynamic, Allocator>>::device_obj(const host_obj& array)
            : length(array.getLength()), capacity(array.getCapacity()) {
        d_data = alloc.allocate(capacity);
        array.toDevice(*this);
    }
    /**
     * Do not launch kernel before copy is finished, so memcpy is async.
     */
    template<class T, class Allocator>
    device_obj<Array<T, Dynamic, Allocator>>::device_obj(const device_obj<Array<T, Dynamic, Allocator>>& obj)
            : length(obj.getLength()), capacity(obj.getCapacity()), alloc(obj.alloc) {
        d_data = alloc.allocate(capacity);
        auto& ctx = Core::CUDAContext::getInstance();
        if constexpr (isTrivial)
            cudaMemcpyAsync(d_data, obj.d_data, length * sizeof(T), cudaMemcpyKind::cudaMemcpyDeviceToDevice, ctx);
        else {
            Array<ElemType, Dynamic> buffer(length);
            cudaMemcpyAsync(buffer.data(), obj.d_data, length * sizeof(ElemType), cudaMemcpyKind::cudaMemcpyDeviceToHost, ctx);
            ctx.wait();
            Array<ElemType, Dynamic> buffer1 = buffer;
            cudaMemcpyAsync(d_data, buffer1.data(), length * sizeof(ElemType), cudaMemcpyKind::cudaMemcpyHostToDevice, ctx);
            buffer.get_allocator().deallocate(buffer.release(), length);
            ctx.wait();
            buffer1.get_allocator().deallocate(buffer1.release(), length);
        }
    }

    template<class T, class Allocator>
    device_obj<Array<T, Dynamic, Allocator>>::device_obj(device_obj<Array<T, Dynamic, Allocator>>&& obj) noexcept
            : d_data(obj.d_data), length(obj.length), capacity(obj.capacity), alloc(std::move(obj.alloc)) {
        obj.d_data = nullptr;
        obj.length = obj.capacity = 0;
    }

    template<class T, class Allocator>
    device_obj<Array<T, Dynamic, Allocator>>::~device_obj() {
        if constexpr (!isTrivial) {
            Array<ElemType, Dynamic> buffer(length);
            auto& ctx = Core::CUDAContext::getInstance();
            cudaMemcpyAsync(buffer.data(), d_data, length * sizeof(ElemType), cudaMemcpyKind::cudaMemcpyDeviceToHost, ctx);
            ctx.wait();
        }
        alloc.deallocate(d_data, capacity);
        d_data = nullptr;
        length = capacity = 0;
    }

    template<class T, class Allocator>
    device_obj<Array<T, Dynamic, Allocator>>&
    device_obj<Array<T, Dynamic, Allocator>>::operator=(device_obj<Array<T, Dynamic, Allocator>> obj) noexcept {
        swap(obj);
        return *this;
    }
    /**
     * Synchronization is expected before this, copy is async to task stream
     */
    template<class T, class Allocator>
    typename device_obj<Array<T, Dynamic, Allocator>>::PlainHostObj
    device_obj<Array<T, Dynamic, Allocator>>::toPlainHost() const {
        PlainHostObj result(getLength());
        check(cudaMemcpy(result.data(), (void*)d_data, getLength() * sizeof(ElemType), cudaMemcpyKind::cudaMemcpyDeviceToHost));
        return result;
    }
    /**
     * Synchronization is expected before this, copy is async to task stream
     */
    template<class T, class Allocator>
    Array<T, Dynamic, Allocator> device_obj<Array<T, Dynamic, Allocator>>::toHost() const {
        host_obj result(length);
        if constexpr (isTrivial)
            check(cudaMemcpy(result.data(), d_data, length * sizeof(ElemType), cudaMemcpyKind::cudaMemcpyDeviceToHost));
        else {
            const auto buffer = toPlainHost();
            for (size_t i = 0; i < length; ++i)
                result[i] = buffer[i].getDerived().toHost();
        }
        return result;
    }

    template<class T, class Allocator>
    inline void device_obj<Array<T, Dynamic, Allocator>>::toHost(host_obj& obj) const {
        obj = toHost();
    }
    /**
     * Synchronization is expected before this, copy is async to task stream
     */
    template<class T, class Allocator>
    void device_obj<Array<T, Dynamic, Allocator>>::reserve(size_t size) {
        assert(size > getCapacity());
        const auto buffer = toPlainHost();
        alloc.deallocate(d_data, capacity);
        d_data = alloc.allocate(size);
        check(cudaMemcpy(d_data, buffer.data(), getLength() * sizeof(ElemType), cudaMemcpyKind::cudaMemcpyHostToDevice));
        capacity = size;
        check(cudaStreamSynchronize(nullptr));
    }
    /**
     * Synchronization is expected before this, copy is async to task stream
     */
    template<class T, class Allocator>
    template<class... Args>
    void device_obj<Array<T, Dynamic, Allocator>>::resize(size_t size, Args&&... args) {
        if (size == length)
            return;
        if (capacity < size)
            reserve(size);

        if (length > size) {
            if constexpr (!isTrivial) {
                const size_t delta = length - size;
                Array<ElemType, Dynamic> buffer{};
                buffer.reserve(delta);
                check(cudaMemcpy(buffer.data(), d_data + size, delta * sizeof(ElemType), cudaMemcpyKind::cudaMemcpyDeviceToHost));
                buffer.setLength(delta);
            }
        }
        else {
            const size_t delta = size - length;
            Array<ElemType, Dynamic> buffer(delta);
            for (size_t i = 0; i < delta; ++i)
                buffer.get_allocator().construct(buffer.data() + i, std::forward<Args>(args)...);
            check(cudaMemcpy(d_data + length, buffer.data(), delta * sizeof(ElemType), cudaMemcpyKind::cudaMemcpyHostToDevice));
            buffer.get_allocator().deallocate(buffer.release(), delta);
            check(cudaStreamSynchronize(nullptr));
        }
        length = size;
    }

    template<class T, class Allocator>
    void device_obj<Array<T, Dynamic, Allocator>>::swap(device_obj& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(d_data, obj.d_data);
        std::swap(length, obj.length);
        std::swap(capacity, obj.capacity);
    }

    template<class T, class Allocator>
    typename device_obj<Array<T, Dynamic, Allocator>>::pointer
    device_obj<Array<T, Dynamic, Allocator>>::release() noexcept {
        auto* copy = d_data;
        d_data = nullptr;
        length = capacity = 0;
        return copy;
    }

    template<class T, class Allocator>
    inline auto Array<T, Dynamic, Allocator>::toDevice() const {
        return device_obj<This>(*this);
    }

    template<class T, class Allocator>
    inline auto Array<T, Dynamic, Allocator>::toDeviceAsync() const {
        device_obj<This> result{};
        result.reserve(getCapacity());
        toDeviceAsync(result);
        return device_obj<This>(std::move(result));
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::toDevice(device_obj<This>& obj) const {
        using ElemType = typename device_obj<This>::ElemType;
        constexpr bool isTrivial = device_obj<This>::isTrivial;
        const size_t length = getLength();
        obj.resize(length);

        auto& ctx = Core::CUDAContext::getInstance();
        if constexpr (isTrivial) {
            cudaMemcpyAsync(obj.data(), this->data(), length * sizeof(ElemType), cudaMemcpyKind::cudaMemcpyHostToDevice, ctx);
            ctx.wait();
        }
        else {
            Array<ElemType, Dynamic> buffer(length);
            for (size_t i = 0; i < length; ++i)
                buffer[i] = this->operator[](i).toDevice();
            cudaMemcpyAsync(obj.data(), buffer.data(), length * sizeof(ElemType), cudaMemcpyKind::cudaMemcpyHostToDevice, ctx);
            ctx.wait();
            buffer.get_allocator().deallocate(buffer.release(), length);
        }
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::toDeviceAsync(device_obj<This>& obj) const {
        static_assert(device_obj<This>::isTrivial, "[Error]: Do not pass non trivial elements to device asynchronously");
        using ElemType = typename device_obj<This>::ElemType;
        const size_t length = getLength();
        obj.resize(length);

        auto& ctx = Core::CUDAContext::getInstance();
        cudaMemcpyAsync(obj.data(), this->data(), length * sizeof(ElemType), cudaMemcpyKind::cudaMemcpyHostToDevice, ctx);
    }
}

namespace Physica {
    template<class T, size_t Length, class Allocator>
    class Traits<Core::device_obj<Core::Array<T, Length, Allocator>>> {
    public:
        using AllocatorType = Core::DeviceAllocator<T>;
        using ElemType = typename AllocatorType::value_type;
        constexpr static size_t SizeAtCompile = Length;
    };
}
