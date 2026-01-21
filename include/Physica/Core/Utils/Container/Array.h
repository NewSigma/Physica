/*
 * Copyright 2020-2026 Weibo He.
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

#include <array>
#include "Physica/CRCoro.h"
#include "Physica/Core/Utils/Allocator/HostAllocator.h"
#include "Physica/Core/Utils/CUDA/device_obj.h"
#include "ArrayImpl/ArrayBase.h"

namespace Physica {
    template<class T, size_t Length = Dynamic, class Allocator = HostAllocator<T>> class Array;
    template<class T> class PageLockedAllocator;
    /**
     * \class Array unifies std::array and std::vector with following features:
     * 1. Customized extensions;
     * 2. Array<bool> is not bitwise;
     */
    template<class T, size_t Length, class Allocator>
    class Array : public ArrayBase<Array<T, Length, Allocator>, Allocator>
                , public CRCoro<Array<T, Length, Allocator>> {
        using This = Array<T, Length, Allocator>;
        using Base = ArrayBase<This, Allocator>;
        using IndexType = Array<size_t, Length>;
        static_assert(std::is_default_constructible<T>::value, "[Error]: Expect default constructible T");
        static_assert(!std::same_as<Allocator, PageLockedAllocator<T>>, "[Error]: Page locked array can not have fixed size");
    private:
        constexpr static size_t Align = std::allocator_traits<Allocator>::Align;

        // We have to use such a verbose alignment since GCC 14.2 complains Align maybe 0.
        alignas(std::max(alignof(T), Align)) std::array<T, Length> arr;
        [[no_unique_address]] Allocator alloc;
    public:
        Array() = default;
        __host__ __device__ explicit Array(size_t length, auto&&... args);
        __host__ __device__ Array(std::initializer_list<T> list);
        template<size_t OtherLength, class OtherAlloc>
        Array(const Array<T, OtherLength, OtherAlloc>& other) noexcept;
        Array(const This&) = default;
        Array(This&&) noexcept = default;
        ~Array() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        template<size_t I>
        [[nodiscard]] constexpr auto&& get(this auto&&) noexcept;

        __host__ __device__ void insert(size_t, auto&&...) = delete;
        __host__ __device__ void resize(size_t length, auto&&... args) noexcept;
        __host__ __device__ void reserve([[maybe_unused]] size_t size) noexcept { assert(size == Length); }
        __host__ __device__ void zeros() noexcept;

        [[nodiscard]] auto toDevice() const;
        [[nodiscard]] auto toDeviceAsync() const;
        void toDevice(device_obj<This>& obj) const;
        void toDeviceAsync(device_obj<This>& obj) const;

        using Base::read;
        using Base::write;
        __host__ __device__ void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ auto* data(this auto&& self) noexcept { return self.arr.data(); }
        [[nodiscard]] __host__ __device__ constexpr static size_t size() noexcept { return Length; }
        [[nodiscard]] __host__ __device__ constexpr static size_t getLength() noexcept { return Length; }
        [[nodiscard]] __host__ __device__ constexpr static size_t getCapacity() noexcept { return Length; }
        [[nodiscard]] auto get_allocator() const noexcept { return alloc; }
        /* Static members */
        [[nodiscard]] __host__ __device__ static This read(size_t length, const T* __restrict p) noexcept;
        [[nodiscard]] static size_t toIndex1D(const IndexType& __restrict shape, const IndexType& __restrict indices) noexcept;
        [[nodiscard]] static IndexType toIndexND(const IndexType& shape, size_t index) noexcept;
        [[nodiscard]] static This generate(std::invocable<size_t> auto fn);
    };

    template<class T, class Allocator>
    class Array<T, Dynamic, Allocator> : public ArrayBase<Array<T, Dynamic, Allocator>, Allocator>
                                             , public CRCoro<Array<T, Dynamic, Allocator>> {
        using This = Array<T, Dynamic, Allocator>;
        using Base = ArrayBase<This, Allocator>;
    private:
        T* arr = nullptr;
        size_t length = 0;
        size_t capacity = 0;
        [[no_unique_address]] Allocator alloc;
    public:
        Array() = default;
        explicit Array(size_t length_, auto&&... args) noexcept(std::is_nothrow_constructible<T, decltype(args)...>::value);
        Array(std::initializer_list<T> list) noexcept;
        template<size_t Length, class OtherAlloc>
        Array(const Array<T, Length, OtherAlloc>& other) noexcept(std::is_nothrow_copy_assignable<T>::value);
        Array(const This& obj) noexcept(std::is_nothrow_copy_constructible<T>::value);
        Array(This&& obj) noexcept;
        ~Array();
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void grow(auto&&... args) noexcept;
        void append(auto&&... args) noexcept;
        void insert(size_t index, auto&&... args) noexcept;
        void resize(size_t size, auto&&... args) noexcept;
        void reserve(size_t size) noexcept;
        void squeeze() noexcept;
        void clear() noexcept;
        [[nodiscard]] T* release() noexcept;
        void zeros() noexcept;
        void read(const T* __restrict p) noexcept;

        [[nodiscard]] auto toDevice() const;
        [[nodiscard]] auto toDeviceAsync() const;
        void toDevice(device_obj<This>& obj) const;
        void toDeviceAsync(device_obj<This>& obj) const;

        using Base::read;
        using Base::write;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard, gnu::returns_nonnull]] __host__ __device__ auto* data(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ size_t size() const noexcept { return length; }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return length; }
        [[nodiscard]] __host__ __device__ size_t getCapacity() const noexcept { return capacity; }
        [[nodiscard]] auto get_allocator() const noexcept { return alloc; }
        /* Static members */
        [[nodiscard]] static This read(size_t length, const T* __restrict p) noexcept;
        [[nodiscard]] static This generate(size_t length, std::invocable<size_t> auto fn);
    private:
        void adjust(size_t size);
        void resizeImpl(size_t size, auto&&... args) noexcept;
        void doubleSpace() noexcept;
    };

    template<class T, size_t Length, class Allocator>
    void swap(Array<T, Length, Allocator>& __restrict array1, Array<T, Length, Allocator>& __restrict array2) noexcept {
        array1.swap(array2);
    }
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
    struct tuple_size<Physica::Array<T, Length, Allocator>> : public integral_constant<std::size_t, Length> {};

    template<std::size_t I, class T, size_t Length, class Allocator>
    struct tuple_element<I, Physica::Array<T, Length, Allocator>> {
        static_assert(Length > 0, "[Error]: Dynamic array is not a tuple");
        using type = T;
    };
}

#include "ArrayImpl/ArrayImpl.h"
