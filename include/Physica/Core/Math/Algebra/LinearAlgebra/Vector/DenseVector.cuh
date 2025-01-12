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

#include "DenseVector.h"
#include "VectorImpl/ContinuousVector.cuh"
#include "Physica/Core/Utils/Container/Array.cuh"

namespace Physica::Core {
    template<Scalar T, size_t Length, class Allocator>
    class device_obj<DenseVector<T, Length, Allocator>>
            : public device_obj<ContinuousVector<DenseVector<T, Length, Allocator>>>
            , public device_obj<Array<T, Length, Allocator>> {
        using Storage = device_obj<Array<T, Length, Allocator>>;
    public:
        using host_obj = DenseVector<T, Length, Allocator>;
        using Base = device_obj<ContinuousVector<host_obj>>;
        using typename Base::ScalarType;
        using Base::SizeAtCompile;
    private:
        using This = device_obj<host_obj>;
    public:
        device_obj() = default;
        using Storage::Storage;
        template<Vector V>
        __host__ __device__ device_obj(const device_obj<V>& v);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        using Base::operator=;
        using Storage::operator=;
        using Storage::operator[];
        /* Operations */
        [[nodiscard]] inline host_obj toHost() const;
        using Base::toHost;
        using Storage::resize;
        using Storage::swap;
        /* Getters */
        using Storage::data;
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Storage::getLength(); }
        [[nodiscard]] __host__ __device__ ScalarType* data_ptr(size_t index) { return data() + index; }
        [[nodiscard]] __host__ __device__ const ScalarType* data_ptr(size_t index) const { return data() + index; }
        /* Static members */
        template<RandomGenerator R>
        [[nodiscard]] inline static This random_uniform(size_t len);
        template<RandomGenerator R>
        [[nodiscard]] inline static This random_normal(size_t len);
        template<class Distribution, RandomGenerator R>
        [[nodiscard]] inline static This random_any(size_t len, Distribution& dist);
    };

    template<Scalar T, size_t Length, class Allocator>
    template<Vector V>
    __host__ __device__ device_obj<DenseVector<T, Length, Allocator>>::device_obj(const device_obj<V>& v)
            : Storage(v.getLength()) {
        v.assign(*this);
    }

    template<Scalar T, size_t Length, class Allocator>
    inline device_obj<DenseVector<T, Length, Allocator>>::host_obj
    device_obj<DenseVector<T, Length, Allocator>>::toHost() const {
        return host_obj(Storage::toHost());
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RandomGenerator R>
    inline device_obj<DenseVector<T, Length, Allocator>> device_obj<DenseVector<T, Length, Allocator>>::random_uniform(
            size_t len) {
        return host_obj::template random_uniform<R>(len).toDevice();
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RandomGenerator R>
    inline device_obj<DenseVector<T, Length, Allocator>> device_obj<DenseVector<T, Length, Allocator>>::random_normal(
            size_t len) {
        return host_obj::template random_normal<R>(len).toDevice();
    }

    template<Scalar T, size_t Length, class Allocator>
    template<class Distribution, RandomGenerator R>
    inline device_obj<DenseVector<T, Length, Allocator>> device_obj<DenseVector<T, Length, Allocator>>::random_any(
            size_t len, Distribution& dist) {
        return host_obj::template random_any<Distribution, R>(len, dist).toDevice();
    }

    template<Scalar T, size_t Length, class Allocator>
    inline device_obj<DenseVector<T, Length, Allocator>> DenseVector<T, Length, Allocator>::toDevice() const {
        return device_obj<DenseVector<T, Length, Allocator>>(*this);
    }
}

namespace Physica {
    template<Scalar T, size_t Length, class Allocator>
    class Traits<Core::device_obj<DenseVector<T, Length, Allocator>>> : public Traits<DenseVector<T, Length, Allocator>> {};
}
