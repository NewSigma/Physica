/*
 * Copyright 2022-2026 Weibo He.
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
#include "VectorImpl/CompactVector.cuh"
#include "Physica/Core/Utils/Container/Array.cuh"

namespace Physica {
    template<Scalar T, size_t Length, class Allocator>
    class device_obj<DenseVector<T, Length, Allocator>>
            : public device_obj<CompactVector<DenseVector<T, Length, Allocator>>>
            , public CRCoro<device_obj<DenseVector<T, Length, Allocator>>> {
        static_assert(!Diffable<T>, "[Error]: Use diffable vector instead");
        using host_obj = DenseVector<T, Length, Allocator>;
        using This = device_obj<host_obj>;
        using Base = device_obj<CompactVector<host_obj>>;
        using Coro = CRCoro<This>;
        using Storage = device_obj<Array<T, Length, Allocator>>;
    public:
        using typename Base::ScalarType;
        using Base::SizeAtCompile;
    private:
        Storage storage;
    public:
        device_obj() = default;
        explicit device_obj(Storage storage) noexcept;
        __host__ __device__ explicit device_obj(size_t length);
        __host__ __device__ device_obj(size_t length, T init);
        __host__ __device__ device_obj(const Vector auto& v);
        device_obj(const host_obj& obj);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        using Base::operator=;
        using Base::operator[];
        /* Operations */
        __host__ __device__ void resize(const Vector auto& x);
        __host__ __device__ void resize(size_t size, auto&&... args) noexcept;
        void reserve(size_t size) noexcept;

        [[nodiscard]] host_obj toHost() const;
        [[nodiscard]] host_obj toHostAsync() const;
        void toHost(host_obj& obj) const;
        void toHostAsync(host_obj& obj) const;
        using Base::toHost;
        using Base::toHostAsync;

        template<RNG R>
        void random_uniform();
        template<RNG R>
        void random_normal();
        template<RNG R>
        void random_any(auto& distribution);

        using Base::read;
        using Base::write;
        void zeros() noexcept;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return storage.getLength(); }
        [[nodiscard]] __host__ __device__ auto* data(this auto&&) noexcept;
        /* Static members */
        template<RNG R>
        [[nodiscard]] static This random_uniform(size_t len);
        template<RNG R>
        [[nodiscard]] static This random_normal(size_t len);
        template<RNG R>
        [[nodiscard]] static This random_any(size_t len, auto& distribution);
    };

    template<Scalar T, size_t Length, class Allocator>
    auto DenseVector<T, Length, Allocator>::toDevice() const {
        return device_obj<DenseVector<T, Length, Allocator>>(*this);
    }

    template<Scalar T, size_t Length, class Allocator>
    auto DenseVector<T, Length, Allocator>::toDeviceAsync() const {
        device_obj<DenseVector<T, Length, Allocator>> result{};
        toDeviceAsync(result);
        return result;
    }
}

namespace Physica {
    template<Scalar T, size_t Length, class Allocator>
    class Traits<device_obj<DenseVector<T, Length, Allocator>>> : public Traits<DenseVector<T, Length, Allocator>> {};
}

#include "VectorImpl/DenseVectorImpl.cuh"
