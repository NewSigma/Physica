/*
 * Copyright 2020-2024 Weibo He.
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

#include <Physica/Core/MultiPrecision/Scalar.h>
#include <Physica/Utils/Container/Array/Array.h>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixOption.h>
#include "VectorImpl/ContinuousVector.h"

namespace Physica::Core {
    template<class T, size_t Length = Dynamic, class Allocator = Utils::HostAllocator<T>>
    class Vector : public ContinuousVector<Vector<T, Length, Allocator>>, public Utils::Array<T, Length, Allocator> {
        using This = Vector<T, Length, Allocator>;
    public:
        using Base = ContinuousVector<This>;
        using Storage = Utils::Array<T, Length, Allocator>;
        using device_obj_type = device_obj<This>;
        using typename Base::ScalarType;
        using typename Base::PlainScalar;
        using typename Base::ColMatrix;
        using typename Base::RowMatrix;
        using Base::SizeAtCompile;
        using Base::isReverseDiff;
    public:
        Vector() = default;
        explicit Vector(size_t length_);
        Vector(size_t length_, const T& value);
        Vector(std::initializer_list<T> list);
        Vector(Storage array) noexcept;
        template<class Derived>
        Vector(const RValueVector<Derived>& v);
        Vector(const This&) = default;
        Vector(This&&) noexcept = default;
        ~Vector() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        using Base::operator=;
        using Storage::operator[];
        /* Operations */
        using Storage::resize;
        using Base::random_uniform;
        using Base::random_normal;
        using Base::random_any;
        [[nodiscard]] Vector copy() const;
        [[nodiscard]] inline device_obj<This> toDevice() const;
        using Base::toDevice;
        using Base::toDeviceAsync;
        using Storage::swap;
        /* Getters */
        using Storage::getLength;
        using Storage::data;
        [[nodiscard]] __host__ __device__ ScalarType* data_ptr(size_t index) { return data() + index; }
        [[nodiscard]] __host__ __device__ const ScalarType* data_ptr(size_t index) const { return data() + index; }
        /* Static members */
        [[nodiscard]] static Vector zeros(size_t len);
        template<class RandomGenerator>
        [[nodiscard]] static Vector random_uniform(size_t len, RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] static Vector random_uniform(const Vector& v1, const Vector& v2, RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] static Vector random_normal(size_t len, RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        [[nodiscard]] static Vector random_any(size_t len, Distribution& dist, RandomGenerator& gen);
        [[nodiscard]] static Vector linspace(T from, T to, size_t count);
    };
}

namespace Physica {
    using namespace Core;

    template<class T, size_t Length, class Allocator>
    class Traits<Vector<T, Length, Allocator>> {
    public:
        using ScalarType = T;
        constexpr static size_t SizeAtCompile = Length;

        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = true;
    };
}

namespace std {
    template<class T, size_t Length>
    inline void swap(Physica::Core::Vector<T, Length>& __restrict v1, Physica::Core::Vector<T, Length>& __restrict v2) noexcept {
        v1.swap(v2);
    }
}

#include "VectorImpl/VectorImpl.h"
