/*
 * Copyright 2022-2025 Weibo He.
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

#include "Physica/PlainStruct.h"
#include "Physica/Core/Exception/CUDA/CUDA.cuh"
#include "Physica/Core/Exception/NoImplException.h"
#include "Physica/Core/Utils/Allocator/DeviceAllocator.cuh"
#include "Array.h"

namespace Physica {
    template<class T, size_t Length, class Allocator>
    class device_obj<Array<T, Length, Allocator>> : public Array<T, Length, Allocator> {
        static_assert(Length != Dynamic, "[Error]: Dynamic length is not implemented");
        static_assert(std::is_trivially_copyable<T>::value, "[Error]: Fixed size array with non-trivial element is not supported, it is seldom used on cuda");
        using host_obj = Array<T, Length, Allocator>;
        using This = device_obj<host_obj>;
        using Base = host_obj;
    public:
        using host_obj::host_obj;
        device_obj() = default;
        __host__ __device__ device_obj(const host_obj& obj) : host_obj(obj) {}
        device_obj(const This& obj) = default;
        device_obj(This&& obj) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        __host__ __device__ This& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        /* Operations */
        using Base::resize;
        [[nodiscard]] __host__ __device__ host_obj toHost() const { return *this; }
        __host__ __device__ void toHost(host_obj& obj) const noexcept { obj = *this; }
        __host__ __device__ void toHostAsync(host_obj& obj) const noexcept { toHost(obj); }
        __host__ __device__ void swap(This& __restrict obj) noexcept { host_obj::swap(obj); }
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
        constexpr static bool isTrivial = std::is_trivially_copyable<T>::value;
        using typename Base::allocator_type;
        using typename Base::value_type;
        using typename Base::pointer;
        using typename Base::const_pointer;
        using typename Base::lvalue_reference;
        using typename Base::const_lvalue_reference;
        using typename Base::rvalue_reference;

        using typename Base::ElemType;
        using PlainElemType = std::conditional<isTrivial, ElemType, Physica::PlainStruct<ElemType>>::type;
        using PlainElemAllocator = Allocator::template rebind_alloc<PlainElemType>;
        using PlainHostObj = Array<PlainElemType, Dynamic, PlainElemAllocator>;
    protected:
        using typename Base::FIteType;
        using typename Base::RIteType;
        using typename Base::CFIteType;
        using typename Base::CRIteType;
    private:
        pointer d_data = nullptr;
        size_t length = 0;
        size_t capacity = 0;
        [[no_unique_address]] allocator_type alloc;
    public:
        device_obj() = default;
        template<class... Args>
        __host__ __device__ explicit device_obj(size_t length_, Args&&... args);
        explicit device_obj(const host_obj& array);
        device_obj(const This& obj);
        device_obj(This&& obj) noexcept;
        __host__ __device__ ~device_obj();
        /* Operators */
        This& operator=(device_obj obj) noexcept { swap(obj); return *this; }
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
        [[nodiscard]] host_obj toHostAsync() const;
        inline void toHost(host_obj& obj) const;
        void toHostAsync(host_obj& obj) const;

        void reserve(size_t size);
        template<class... Args> void resize(size_t size, Args&&... args);
        void swap(This& __restrict obj) noexcept;
        [[nodiscard]] pointer release() noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t size() const noexcept { return getLength(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return length; }
        [[nodiscard]] __host__ __device__ size_t getCapacity() const noexcept { return capacity; }
        using Base::empty;
        [[nodiscard]] __host__ __device__ pointer data() noexcept { return d_data; }
        [[nodiscard]] __host__ __device__ const_pointer data() const noexcept { return d_data; }
    private:
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
    template<class... Args>
    __host__ __device__ device_obj<Array<T, Dynamic, Allocator>>::device_obj(size_t length_, Args&&... args) {
        if constexpr (IsHost())
            host_obj(length_, std::forward<Args>(args)...).toDevice(*this);
        else
            unreachable(); // Marked __device__ to silence warnings
    }

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
        auto& ctx = CUDAContext::getInstance();
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
    __host__ __device__ device_obj<Array<T, Dynamic, Allocator>>::~device_obj() {
        if constexpr (IsHost()) {
            if constexpr (!isTrivial) {
                if (length != 0) {
                    Array<ElemType, Dynamic> buffer(length);
                    auto& ctx = CUDAContext::getInstance();
                    cudaMemcpyAsync(buffer.data(), d_data, length * sizeof(ElemType), cudaMemcpyKind::cudaMemcpyDeviceToHost, ctx);
                    ctx.wait();
                }
            }
            alloc.deallocate(d_data, capacity);
            d_data = nullptr;
            length = capacity = 0;
        }
        else
            unreachable(); // Marked __device__ to silence warnings
    }

    template<class T, class Allocator>
    auto device_obj<Array<T, Dynamic, Allocator>>::toPlainHost() const -> PlainHostObj {
        assert(getLength() > 0 && "[Error]: Unnecessary memcpy a empty array");
        PlainHostObj result(getLength());
        check(cudaMemcpyAsync(result.data(), (void*)d_data, getLength() * sizeof(ElemType), cudaMemcpyKind::cudaMemcpyDeviceToHost, CUDAContext::getInstance()));
        return result;
    }

    template<class T, class Allocator>
    Array<T, Dynamic, Allocator> device_obj<Array<T, Dynamic, Allocator>>::toHost() const {
        host_obj result(length);
        toHost(result);
        return result;
    }

    template<class T, class Allocator>
    inline void device_obj<Array<T, Dynamic, Allocator>>::toHost(host_obj& obj) const {
        toHostAsync(obj);
        CUDAContext::getInstance().wait();
    }

    template<class T, class Allocator>
    auto device_obj<Array<T, Dynamic, Allocator>>::toHostAsync() const -> host_obj {
        host_obj result(length);
        toHostAsync(result);
        return result;
    }

    template<class T, class Allocator>
    void device_obj<Array<T, Dynamic, Allocator>>::toHostAsync(host_obj& obj) const {
        obj.reserve(length);
        if constexpr (isTrivial)
            check(cudaMemcpyAsync(obj.data(), d_data, length * sizeof(ElemType), cudaMemcpyKind::cudaMemcpyDeviceToHost, CUDAContext::getInstance()));
        else {
            const auto buffer = toPlainHost();
            for (size_t i = 0; i < length; ++i)
                buffer[i].getDerived().toHostAsync(obj[i]);
        }
        obj.setLength(length);
    }

    template<class T, class Allocator>
    void device_obj<Array<T, Dynamic, Allocator>>::reserve(size_t size) {
        assert(size > getCapacity());
        if (getCapacity() == 0 || getLength() == 0)
            d_data = alloc.allocate(size);
        else {
            const auto buffer = toPlainHost();
            alloc.deallocate(d_data, capacity);
            d_data = alloc.allocate(size);
            cudaMemcpyAsync(d_data, buffer.data(), getLength() * sizeof(ElemType), cudaMemcpyKind::cudaMemcpyHostToDevice, CUDAContext::getInstance());
        }
        capacity = size;
    }

    template<class T, class Allocator>
    template<class... Args>
    void device_obj<Array<T, Dynamic, Allocator>>::resize(size_t size, Args&&... args) {
        if (size == length)
            return;
        if (capacity < size)
            reserve(size);

        auto& ctx = CUDAContext::getInstance();
        if (length > size) {
            if constexpr (!isTrivial) {
                const size_t delta = length - size;
                Array<ElemType, Dynamic> buffer{};
                buffer.reserve(delta);
                check(cudaMemcpyAsync(buffer.data(), d_data + size, delta * sizeof(ElemType), cudaMemcpyKind::cudaMemcpyDeviceToHost, ctx));
                buffer.setLength(delta);
            }
        }
        else if constexpr (!Base::template isTrivialDefaultConstruct<Args...>()) {
            const size_t delta = size - length;
            Array<ElemType, Dynamic> buffer(delta);
            for (size_t i = 0; i < delta; ++i)
                buffer.get_allocator().construct(buffer.data() + i, std::forward<Args>(args)...);
            check(cudaMemcpyAsync(d_data + length, buffer.data(), delta * sizeof(ElemType), cudaMemcpyKind::cudaMemcpyHostToDevice, ctx));
            buffer.get_allocator().deallocate(buffer.release(), delta);
        }
        length = size;
    }

    template<class T, class Allocator>
    void device_obj<Array<T, Dynamic, Allocator>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(d_data, obj.d_data);
        std::swap(length, obj.length);
        std::swap(capacity, obj.capacity);
    }

    template<class T, class Allocator>
    auto device_obj<Array<T, Dynamic, Allocator>>::release() noexcept -> pointer {
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
        using ElemType = device_obj<This>::ElemType;
        constexpr bool isTrivial = device_obj<This>::isTrivial;
        const size_t length = getLength();
        obj.resize(length);

        auto& ctx = CUDAContext::getInstance();
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
        using ElemType = device_obj<This>::ElemType;
        const size_t length = getLength();
        obj.resize(length);

        auto& ctx = CUDAContext::getInstance();
        cudaMemcpyAsync(obj.data(), this->data(), length * sizeof(ElemType), cudaMemcpyKind::cudaMemcpyHostToDevice, ctx);
    }
}

namespace Physica {
    template<class T, size_t Length, class Allocator>
    class Traits<device_obj<Array<T, Length, Allocator>>> {
    public:
        using AllocatorType = DeviceAllocator<T>;
        using ElemType = AllocatorType::value_type;
        constexpr static size_t SizeAtCompile = Length;
    };
}
