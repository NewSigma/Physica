/*
 * Copyright 2022-2024 WeiBo He.
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

namespace Physica::Core {
    namespace Internal {
        template<class T, size_t Length, size_t MaxLength>
        class Traits<Core::device_obj<Vector<T, Length, MaxLength>>> : public Traits<Vector<T, Length, MaxLength>> {};
    }

    template<class T, size_t Length, size_t MaxLength>
    class device_obj<Vector<T, Length, MaxLength>>
            : public device_obj<ContinuousVector<Vector<T, Length, MaxLength>>>
            , public Utils::device_obj<Utils::Array<T, Length, MaxLength>> {
        using Storage = Utils::device_obj<Utils::Array<T, Length, MaxLength>>;
    public:
        using host_obj = Vector<T, Length, MaxLength>;
        using Base = device_obj<ContinuousVector<host_obj>>;
        using typename Base::ScalarType;
    public:
        device_obj() = default;
        using Storage::Storage;
        template<class Derived>
        device_obj(const device_obj<RValueVector<Derived>>& v);
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj<Vector<T, Length, MaxLength>>& operator=(device_obj<Vector<T, Length, MaxLength>> obj) noexcept;
        using Base::operator=;
        using Storage::operator=;
        using Storage::operator[];
        /* Operations */
        [[nodiscard]] inline host_obj toHost() const;
        using Base::toHost;
        using Storage::resize;
        using Storage::swap;
        /* Getters */
        using Storage::getLength;
        using Storage::data;
        [[nodiscard]] __host__ __device__ ScalarType* data_ptr(size_t index) { return data() + index; }
        [[nodiscard]] __host__ __device__ const ScalarType* data_ptr(size_t index) const { return data() + index; }
    };

    template<class T, size_t Length, size_t MaxLength>
    template<class Derived>
    device_obj<Vector<T, Length, MaxLength>>::device_obj(const device_obj<RValueVector<Derived>>& v) : Storage(v.getLength()) {
        v.getDerived().assignTo(*this);
    }

    template<class T, size_t Length, size_t MaxLength>
    device_obj<Vector<T, Length, MaxLength>>&
    device_obj<Vector<T, Length, MaxLength>>::operator=(device_obj<Vector<T, Length, MaxLength>> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class T, size_t Length, size_t MaxLength>
    inline typename device_obj<Vector<T, Length, MaxLength>>::host_obj
    device_obj<Vector<T, Length, MaxLength>>::toHost() const {
        return host_obj(Storage::toHost());
    }

    template<class T, size_t Length, size_t MaxLength, class Allocator>
    inline device_obj<Vector<T, Length, MaxLength, Allocator>> Vector<T, Length, MaxLength, Allocator>::toDevice() const {
        return device_obj<Vector<T, Length, MaxLength>>(*this);
    }
}
