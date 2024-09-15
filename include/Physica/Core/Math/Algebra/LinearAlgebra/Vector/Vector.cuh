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

#include "Vector.h"
#include "VectorImpl/ContinuousVector.cuh"
#include <Physica/Core/Utils/Container/Array.cuh>

namespace Physica::Core {
    template<class T, size_t Length>
    class device_obj<Vector<T, Length>>
            : public device_obj<ContinuousVector<Vector<T, Length>>>
            , public device_obj<Array<T, Length>> {
        using Storage = device_obj<Array<T, Length>>;
    public:
        using host_obj = Vector<T, Length>;
        using Base = device_obj<ContinuousVector<host_obj>>;
        using typename Base::ScalarType;
        using Base::SizeAtCompile;
    private:
        using This = device_obj<host_obj>;
    public:
        device_obj() = default;
        using Storage::Storage;
        template<class Derived>
        __host__ __device__ device_obj(const device_obj<RValueVector<Derived>>& v);
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
        template<class RandomGenerator>
        [[nodiscard]] inline static This random_uniform(size_t len, RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] inline static This random_normal(size_t len, RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        [[nodiscard]] inline static This random_any(size_t len, Distribution& dist, RandomGenerator& gen);
    };

    template<class T, size_t Length>
    template<class Derived>
    __host__ __device__ device_obj<Vector<T, Length>>::device_obj(const device_obj<RValueVector<Derived>>& v) : Storage(v.getLength()) {
        v.getDerived().assignTo(*this);
    }

    template<class T, size_t Length>
    inline typename device_obj<Vector<T, Length>>::host_obj
    device_obj<Vector<T, Length>>::toHost() const {
        return host_obj(Storage::toHost());
    }

    template<class T, size_t Length>
    template<class RandomGenerator>
    inline device_obj<Vector<T, Length>> device_obj<Vector<T, Length>>::random_uniform(
            size_t len, RandomGenerator& gen) {
        return host_obj::random_uniform(len, gen).toDevice();
    }

    template<class T, size_t Length>
    template<class RandomGenerator>
    inline device_obj<Vector<T, Length>> device_obj<Vector<T, Length>>::random_normal(
            size_t len, RandomGenerator& gen) {
        return host_obj::random_normal(len, gen).toDevice();
    }

    template<class T, size_t Length>
    template<class Distribution, class RandomGenerator>
    inline device_obj<Vector<T, Length>> device_obj<Vector<T, Length>>::random_any(
            size_t len, Distribution& dist, RandomGenerator& gen) {
        return host_obj::random_any(len, dist, gen).toDevice();
    }

    template<class T, size_t Length, class Allocator>
    inline device_obj<Vector<T, Length, Allocator>> Vector<T, Length, Allocator>::toDevice() const {
        return device_obj<Vector<T, Length>>(*this);
    }
}

namespace Physica {
    using namespace Core;

    template<class T, size_t Length>
    class Traits<Core::device_obj<Vector<T, Length>>> : public Traits<Vector<T, Length>> {};
}
