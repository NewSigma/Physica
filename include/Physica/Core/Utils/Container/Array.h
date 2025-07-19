/*
 * Copyright 2020-2025 Weibo He.
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

#include <cassert>
#include <cstdlib>
#include <initializer_list>
#include "Physica/CRCoro.h"
#include "Physica/Core/Utils/Allocator/HostAllocator.h"
#include "Physica/Core/Utils/CUDA/device_obj.cuh"
#include "ArrayImpl/ArrayBase.h"

namespace Physica {
    template<class T, size_t Length = Dynamic, class Allocator = HostAllocator<T>> class Array;
    template<class T> class PageLockedAllocator;

    template<class T, size_t Length, class Allocator>
    class Array : public ArrayBase<Array<T, Length, Allocator>, Allocator>
                , public CRCoro<Array<T, Length, Allocator>> {
        static_assert(!std::is_same_v<Allocator, PageLockedAllocator<T>>, "[Error]: Page locked array can not have fixed size");
        using This = Array<T, Length, Allocator>;
    public:
        using Base = ArrayBase<This, Allocator>;
        using typename Base::ElemType;
        using typename Base::allocator_type;
        using typename Base::pointer;
        using typename Base::const_pointer;
        using typename Base::lvalue_reference;
        using typename Base::const_lvalue_reference;
        using typename Base::rvalue_reference;
    private:
        constexpr static size_t Align = std::allocator_traits<Allocator>::Align;

    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wattributes"
        alignas(Align) T arr[Length];
    #pragma GCC diagnostic pop
        [[no_unique_address]] allocator_type alloc;
    public:
        Array() = default;
        __host__ __device__ explicit Array(size_t length, auto&&... args);
        __host__ __device__ Array(std::initializer_list<T> list);
        template<size_t OtherLength, class OtherAlloc>
        Array(const Array<T, OtherLength, OtherAlloc>& other);
        Array(const This&) = default;
        Array(This&&) noexcept = default;
        ~Array() = default;
        /* Operators */
        __host__ __device__ This& operator=(This array) noexcept { swap(array); return *this; }
        /* Operations */
        __host__ __device__ void insert(const T&, size_t) { assert(false); }
        __host__ __device__ void reserve([[maybe_unused]] size_t size) { assert(size == Length); }
        __host__ __device__ void resize(size_t length, auto&&... args) noexcept;
        __host__ __device__ void zeros() noexcept;

        [[nodiscard]] auto toDevice() const;
        [[nodiscard]] auto toDeviceAsync() const;
        void toDevice(device_obj<This>& obj) const;
        void toDeviceAsync(device_obj<This>& obj) const;

        using Base::read;
        using Base::write;
        __host__ __device__ void swap(Array& __restrict array) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr static size_t size() { return Length; }
        [[nodiscard]] __host__ __device__ constexpr static size_t getLength() { return Length; }
        [[nodiscard]] __host__ __device__ constexpr static size_t getCapacity() { return Length; }
        [[nodiscard]] __host__ __device__ pointer data() noexcept { return arr; }
        [[nodiscard]] __host__ __device__ const_pointer data() const noexcept { return arr; }
        [[nodiscard]] allocator_type get_allocator() const noexcept { return alloc; }
        /* Setters */
        void setLength([[maybe_unused]] size_t size) { assert(size == Length); }
        /* Static members */
        [[nodiscard]] __host__ __device__ static This read(size_t length, const T* __restrict p);
    };

    template<class T, class Allocator>
    class Array<T, Dynamic, Allocator> : public ArrayBase<Array<T, Dynamic, Allocator>, Allocator>
                                       , public CRCoro<Array<T, Dynamic, Allocator>> {
        using This = Array<T, Dynamic, Allocator>;
    public:
        using Base = ArrayBase<This, Allocator>;
        using typename Base::ElemType;
        using typename Base::allocator_type;
        using typename Base::pointer;
        using typename Base::const_pointer;
        using typename Base::lvalue_reference;
        using typename Base::const_lvalue_reference;
        using typename Base::rvalue_reference;
    private:
        using Base::getDerived;
        constexpr static size_t MinDeltaSpace = 1024;
        constexpr static size_t Align = std::allocator_traits<Allocator>::Align;
    private:
        pointer arr = nullptr;
        size_t length = 0;
        size_t capacity = 0;
        [[no_unique_address]] allocator_type alloc;
    public:
        Array() = default;
        explicit Array(size_t length_, auto&&... args);
        Array(std::initializer_list<T> list);
        template<size_t Length, class OtherAlloc>
        Array(const Array<T, Length, OtherAlloc>& other);
        Array(const This& obj);
        Array(This&& array) noexcept;
        ~Array();
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void grow(auto&&... args);
        void append(auto&&... args);
        void insert(size_t index, auto&&... args);
        void reserve(size_t size);
        void resize(size_t size, auto&&... args);
        void squeeze();
        void increase(size_t size);
        void decrease(size_t size);
        void clear() noexcept;
        [[nodiscard]] pointer release() noexcept;
        void doubleSpace();
        void zeros() noexcept;

        [[nodiscard]] auto toDevice() const;
        [[nodiscard]] auto toDeviceAsync() const;
        void toDevice(device_obj<This>& obj) const;
        void toDeviceAsync(device_obj<This>& obj) const;

        using Base::read;
        using Base::write;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ pointer data() noexcept;
        [[nodiscard]] __host__ __device__ const_pointer data() const noexcept;
        [[nodiscard]] __host__ __device__ size_t size() const noexcept { return length; }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return length; }
        [[nodiscard]] __host__ __device__ size_t getCapacity() const noexcept { return capacity; }
        [[nodiscard]] allocator_type get_allocator() const noexcept { return alloc; }
        /* Setters */
        void setLength(size_t size);
        /* Static members */
        [[nodiscard]] static This read(size_t length, const T* __restrict p);
    private:
        void resizeImpl(size_t size, auto&&... args);
    };
}

namespace Physica {
    template<class T, size_t Length, class Allocator>
    class Traits<Array<T, Length, Allocator>> {
    public:
        using ElemType = T;
        constexpr static size_t SizeAtCompile = Length;
        using AllocatorType = Allocator;
    };
}

namespace std {
    template<class T, size_t Length, class Allocator>
    void swap(Physica::Array<T, Length, Allocator>& __restrict array1,
                     Physica::Array<T, Length, Allocator>& __restrict array2) noexcept {
        array1.swap(array2);
    }
}

#include "ArrayImpl/ArrayImpl.h"
