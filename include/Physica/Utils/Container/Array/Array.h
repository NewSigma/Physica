/*
 * Copyright 2020-2022 WeiBo He.
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
#include "DynamicArrayBase.h"
#include "Physica/Utils/Container/HostAllocator.h"
#include "Physica/Utils/CUDA/device_obj.cuh"

namespace Physica::Utils {
    constexpr size_t Dynamic = 0;
    /**
     * Linear storage container.
     * This class is designed to be a dynamic or fixed array, in detail:
     * 1.A fixed length array.(same to std::array)
     * 2.A dynamic array but its max size is fixed.
     * 3.A dynamic array whose length and max size are dynamic.(same to std::vector)
     * Arrays can be casted to each other.
     *
     * Note:
     * If \T is a complex type, \T must have its copy and move constructors defined.
     *
     * Optimize:
     * 1. Use the end pointer of arr instead of length may improve performance.
     */
    template<class T, size_t Length = Dynamic, size_t Capacity = Length, class Allocator = HostAllocator<T>>
    class Array;

    namespace Internal {
        template<class T, size_t Length, size_t Capacity, class Allocator>
        class Traits<Array<T, Length, Capacity, Allocator>> {
        public:
            using ValueType = T;
            constexpr static size_t ArrayLength = Length;
            constexpr static size_t ArrayCapacity = Capacity;
            using AllocatorType = Allocator;
        };
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    class Array : public Internal::ArrayBase<Array<T, Length, Capacity, Allocator>, Allocator> {
        static_assert(Length == Capacity, "[Error]: Capacity of fixed array must equals to Length.");
        static_assert(sizeof(T) * Length <= (1U << 16U), "[Warning]: Allocate large fixed array on stack is not recommanded");
    public:
        using This = Array<T, Length, Capacity, Allocator>;
        using Base = Internal::ArrayBase<This, Allocator>;
        using typename Base::allocator_type;
        using typename Base::ValueType;
        using typename Base::pointer;
        using typename Base::const_pointer;
        using typename Base::lvalue_reference;
        using typename Base::const_lvalue_reference;
        using typename Base::rvalue_reference;
        constexpr static size_t ArrayLength = Length;
        constexpr static size_t ArrayCapacity = Capacity;
    private:
        T arr[Length];
        allocator_type alloc;
    public:
        Array() = default;
        template<class... Args>
        explicit Array(size_t length_, Args... args);
        Array(std::initializer_list<T> list);
        Array(const Array& array);
        Array(Array&& array) noexcept;
        ~Array() = default;
        /* Operators */
        __host__ __device__ Array& operator=(Array array) noexcept { swap(array); return *this; }
        /* Operations */
        Array<T, Dynamic, Dynamic, Allocator> subArray(size_t from, size_t to);
        Array<T, Dynamic, Dynamic, Allocator> subArray(size_t from) { return subArray(from, Length); }
        Array<T, Dynamic, Dynamic, Allocator> cut(size_t from);
        void insert(const T&, size_t) { assert(false); }
        void reserve([[maybe_unused]] size_t size) { assert(size == Capacity); }
        template<class... Args>
        void resize([[maybe_unused]] size_t size, [[maybe_unused]] Args... args) { assert(size == Length); }
        __host__ __device__ void swap(Array& array) noexcept;
        [[nodiscard]] inline device_obj<This> toDevice() const;
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr static size_t size() { return Length; }
        [[nodiscard]] __host__ __device__ constexpr static size_t getLength() { return Length; }
        [[nodiscard]] __host__ __device__ constexpr static size_t getCapacity() { return Capacity; }
        [[nodiscard]] pointer data() noexcept { return arr; }
        [[nodiscard]] const_pointer data() const noexcept { return arr; }
        [[nodiscard]] allocator_type get_allocator() const noexcept { return alloc; }
        /* Setters */
        void setLength([[maybe_unused]] size_t size) { assert(size == Length); }
    };

    template<class T, size_t Capacity, class Allocator>
    class Array<T, Dynamic, Capacity, Allocator>
        : public Internal::DynamicArrayBase<Array<T, Dynamic, Capacity, Allocator>, Allocator> {
    public:
        using Base = Internal::DynamicArrayBase<Array<T, Dynamic, Capacity, Allocator>, Allocator>;
        using typename Base::ValueType;
        using typename Base::pointer;
        using typename Base::lvalue_reference;
        using typename Base::const_lvalue_reference;
        using typename Base::rvalue_reference;
        constexpr static size_t ArrayLength = Dynamic;
        constexpr static size_t ArrayCapacity = Capacity;
    private:
        using Base::length;
        using Base::arr;
        using Base::alloc;
        using Base::getDerived;
    public:
        Array();
        template<class... Args>
        explicit Array(size_t length_, Args... args);
        Array(std::initializer_list<T> list);
        Array(const Array& array);
        Array(Array&& array) noexcept;
        template<size_t OtherLength, size_t OtherCapacity>
        explicit Array(const Array<T, OtherLength, OtherCapacity, Allocator>& array);
        template<size_t OtherLength, size_t OtherCapacity>
        explicit Array(Array<T, OtherLength, OtherCapacity, Allocator>&& array) noexcept;
        template<class OtherT, size_t OtherLength, size_t OtherCapacity>
        explicit Array(const Array<OtherT, OtherLength, OtherCapacity, Allocator>& array);
        template<class OtherT, size_t OtherLength, size_t OtherCapacity>
        explicit Array(Array<OtherT, OtherLength, OtherCapacity, Allocator>&& array) noexcept;
        ~Array() = default;
        /* Operators */
        Array& operator=(Array array) noexcept { swap(array); return *this; }
        /* Operations */
        Array<T, Dynamic, Dynamic, Allocator> subArray(size_t from, size_t to);
        Array<T, Dynamic, Dynamic, Allocator> subArray(size_t from) { return subArray(from, length); }
        Array<T, Dynamic, Dynamic, Allocator> cut(size_t from);
        inline void append(const_lvalue_reference t);
        inline void append(rvalue_reference t);
        void append(const Array& t);
        void append(Array&& t);
        void reserve(size_t size);
        void swap(Array& array) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t size() const noexcept { return length; }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return length; }
        [[nodiscard]] __host__ __device__ constexpr static size_t getCapacity() { return Capacity; }
    };

    template<class T, class Allocator>
    class Array<T, Dynamic, Dynamic, Allocator>
        : public Internal::DynamicArrayBase<Array<T, Dynamic, Dynamic, Allocator>, Allocator> {
    public:
        using This = Array<T, Dynamic, Dynamic, Allocator>;
        using Base = Internal::DynamicArrayBase<This, Allocator>;
        using typename Base::ValueType;
        using typename Base::pointer;
        using typename Base::lvalue_reference;
        using typename Base::const_lvalue_reference;
        using typename Base::rvalue_reference;
        constexpr static size_t ArrayLength = Dynamic;
        constexpr static size_t ArrayCapacity = Dynamic;
    private:
        using Base::length;
        using Base::arr;
        using Base::alloc;
        using Base::getDerived;
        constexpr static size_t MinDeltaSpace = 1024;
    protected:
        size_t capacity;
    public:
        Array();
        template<class... Args>
        explicit Array(size_t length_, Args... args);
        Array(std::initializer_list<T> list);
        Array(const Array& array);
        Array(Array&& array) noexcept;
        template<size_t OtherLength, size_t OtherCapacity>
        explicit Array(const Array<T, OtherLength, OtherCapacity, Allocator>& array);
        template<size_t OtherLength, size_t OtherCapacity>
        explicit Array(Array<T, OtherLength, OtherCapacity, Allocator>&& array) noexcept;
        template<class OtherT, size_t OtherLength, size_t OtherCapacity>
        explicit Array(const Array<OtherT, OtherLength, OtherCapacity, Allocator>& array);
        template<class OtherT, size_t OtherLength, size_t OtherCapacity>
        explicit Array(Array<OtherT, OtherLength, OtherCapacity, Allocator>&& array) noexcept;
        ~Array() = default;
        /* Operators */
        Array& operator=(Array array) noexcept { swap(array); return *this; }
        /* Operations */
        void append(const_lvalue_reference t);
        void append(rvalue_reference t);
        void append(const Array& t);
        void append(Array&& t);
        template<class... Args>
        void insert(size_t index, Args... args);
        void reserve(size_t size);
        template<class... Args>
        void resize(size_t size, Args... args);
        void squeeze();
        void increase(size_t size);
        void decrease(size_t size);
        void swap(Array& array) noexcept;
        [[nodiscard]] inline device_obj<This> toDevice() const;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t size() const noexcept { return length; }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return length; }
        [[nodiscard]] __host__ __device__ size_t getCapacity() const noexcept { return capacity; }
    private:
        void doubleSpace() { increase(capacity * 2 + (MinDeltaSpace + sizeof(T) - 1) / sizeof(T)); }
    };

    template<class T, size_t Length, size_t Capacity, class Allocator>
    inline void swap(Physica::Utils::Array<T, Length, Capacity, Allocator>& array1,
                     Physica::Utils::Array<T, Length, Capacity, Allocator>& array2) noexcept {
        array1.swap(array2);
    }
}

#include "ArrayImpl.h"

#ifdef PHYSICA_CUDA
    #include "Array.cuh"
#endif
