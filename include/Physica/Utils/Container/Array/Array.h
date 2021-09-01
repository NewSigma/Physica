/*
 * Copyright 2020-2021 WeiBo He.
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
     * Unfinished:
     * 1.Copy, move constructors and assign operators maybe able to accept different specializations.
     * 2.Rename 
     *
     * Optimize:
     * 1. Use more effective allocate strategy to avoid reallocate.
     * 2. Use the end pointer of arr instead of length may improve performance.
     * 
     * Idea:
     * May be extends std::array and std::vector.
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
    class Array : public Internal::ArrayStorage<Array<T, Length, Capacity, Allocator>> {
        static_assert(Length == Capacity, "Capacity of fixed array must equals to Length.");
    public:
        using Base = Internal::ArrayStorage<Array<T, Length, Capacity, Allocator>>;
        using typename Base::ValueType;
        using typename Base::PointerType;
        using typename Base::LValueReferenceType;
        using typename Base::ConstLValueReferenceType;
        using typename Base::RValueReferenceType;
        constexpr static size_t ArrayLength = Length;
        constexpr static size_t ArrayCapacity = Capacity;
    private:
        using Base::arr;
        using Base::alloc;
        using Base::getDerived;
    public:
        Array();
        explicit Array(size_t length_, ConstLValueReferenceType t = T());
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
        Array& operator=(const Array& array);
        Array& operator=(Array&& array) noexcept;
        template<size_t OtherLength, size_t OtherCapacity>
        Array& operator=(const Array<T, OtherLength, OtherCapacity, Allocator>& array);
        template<size_t OtherLength, size_t OtherCapacity>
        Array& operator=(Array<T, OtherLength, OtherCapacity, Allocator>&& array) noexcept;
        template<class OtherT, size_t OtherLength, size_t OtherCapacity>
        Array& operator=(const Array<OtherT, OtherLength, OtherCapacity, Allocator>& array);
        template<class OtherT, size_t OtherLength, size_t OtherCapacity>
        Array& operator=(Array<OtherT, OtherLength, OtherCapacity, Allocator>&& array) noexcept;
        /* Helpers */
        Array<T, Dynamic, Dynamic, Allocator> subArray(size_t from, size_t to);
        Array<T, Dynamic, Dynamic, Allocator> subArray(size_t from) { return subArray(from, Length); }
        Array<T, Dynamic, Dynamic, Allocator> cut(size_t from);
        void insert(const T&, size_t) { assert(false); }
        void reserve([[maybe_unused]] size_t size) { assert(size == Capacity); }
        void resize([[maybe_unused]] size_t size) { assert(size == Length); }
        void resize(size_t size, const T& t);
        void swap(Array& array) noexcept { Base::swap(array); }
        /* Getters */
        [[nodiscard]] constexpr static size_t size() { return Length; }
        [[nodiscard]] constexpr static size_t getLength() { return Length; }
        [[nodiscard]] constexpr static size_t getCapacity() { return Capacity; }
        /* Setters */
        void setLength(size_t size) { assert(size == Length); }
    };

    template<class T, size_t Capacity, class Allocator>
    class Array<T, Dynamic, Capacity, Allocator>
        : public Internal::DynamicArrayBase<Array<T, Dynamic, Capacity, Allocator>> {
    public:
        using Base = Internal::DynamicArrayBase<Array<T, Dynamic, Capacity, Allocator>>;
        using typename Base::ValueType;
        using typename Base::PointerType;
        using typename Base::LValueReferenceType;
        using typename Base::ConstLValueReferenceType;
        using typename Base::RValueReferenceType;
        constexpr static size_t ArrayLength = Dynamic;
        constexpr static size_t ArrayCapacity = Capacity;
    private:
        using Base::length;
        using Base::arr;
        using Base::alloc;
        using Base::getDerived;
    public:
        Array();
        explicit Array(size_t length_, ConstLValueReferenceType t = T());
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
        Array& operator=(const Array& array);
        Array& operator=(Array&& array) noexcept;
        template<size_t OtherLength, size_t OtherCapacity>
        Array& operator=(const Array<T, OtherLength, OtherCapacity, Allocator>& array);
        template<size_t OtherLength, size_t OtherCapacity>
        Array& operator=(Array<T, OtherLength, OtherCapacity, Allocator>&& array) noexcept;
        template<class OtherT, size_t OtherLength, size_t OtherCapacity>
        Array& operator=(const Array<OtherT, OtherLength, OtherCapacity, Allocator>& array);
        template<class OtherT, size_t OtherLength, size_t OtherCapacity>
        Array& operator=(Array<OtherT, OtherLength, OtherCapacity, Allocator>&& array) noexcept;
        /* Helpers */
        Array<T, Dynamic, Dynamic, Allocator> subArray(size_t from, size_t to);
        Array<T, Dynamic, Dynamic, Allocator> subArray(size_t from) { return subArray(from, length); }
        Array<T, Dynamic, Dynamic, Allocator> cut(size_t from);
        inline void append(ConstLValueReferenceType t);
        inline void append(RValueReferenceType t);
        void append(const Array& t);
        void append(Array&& t);
        void reserve(size_t size);
        void swap(Array& array) noexcept;
        /* Getters */
        [[nodiscard]] size_t size() const noexcept { return length; }
        [[nodiscard]] size_t getLength() const noexcept { return length; }
        [[nodiscard]] constexpr static size_t getCapacity() { return Capacity; }
    };

    template<class T, class Allocator>
    class Array<T, Dynamic, Dynamic, Allocator>
        : public Internal::DynamicArrayBase<Array<T, Dynamic, Dynamic, Allocator>> {
    public:
        using Base = Internal::DynamicArrayBase<Array<T, Dynamic, Dynamic, Allocator>>;
        using typename Base::ValueType;
        using typename Base::PointerType;
        using typename Base::LValueReferenceType;
        using typename Base::ConstLValueReferenceType;
        using typename Base::RValueReferenceType;
        constexpr static size_t ArrayLength = Dynamic;
        constexpr static size_t ArrayCapacity = Dynamic;
    private:
        using Base::length;
        using Base::arr;
        using Base::alloc;
        using Base::getDerived;
    protected:
        size_t capacity;
    public:
        Array();
        explicit Array(size_t length_, ConstLValueReferenceType t = T());
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
        Array& operator=(const Array& array);
        Array& operator=(Array&& array) noexcept;
        template<size_t OtherLength, size_t OtherCapacity>
        Array& operator=(const Array<T, OtherLength, OtherCapacity, Allocator>& array);
        template<size_t OtherLength, size_t OtherCapacity>
        Array& operator=(Array<T, OtherLength, OtherCapacity, Allocator>&& array) noexcept;
        template<class OtherT, size_t OtherLength, size_t OtherCapacity>
        Array& operator=(const Array<OtherT, OtherLength, OtherCapacity, Allocator>& array);
        template<class OtherT, size_t OtherLength, size_t OtherCapacity>
        Array& operator=(Array<OtherT, OtherLength, OtherCapacity, Allocator>&& array) noexcept;
        /* Helpers */
        void append(ConstLValueReferenceType t);
        void append(RValueReferenceType t);
        void append(const Array& t);
        void append(Array&& t);
        void reserve(size_t size);
        void resize(size_t size);
        void resize(size_t size, const T& t);
        void squeeze();
        void increase(size_t size);
        void decrease(size_t size);
        void swap(Array& array) noexcept;
        /* Getters */
        [[nodiscard]] size_t size() const noexcept { return length; }
        [[nodiscard]] size_t getLength() const noexcept { return length; }
        [[nodiscard]] size_t getCapacity() const noexcept { return capacity; }
    };

    template<class T, size_t Length, size_t Capacity, class Allocator>
    inline void swap(Array<T, Length, Capacity, Allocator>& array1, Array<T, Length, Capacity, Allocator>& array2) noexcept {
        array1.swap(array2);
    }
}

#include "ArrayImpl.h"
