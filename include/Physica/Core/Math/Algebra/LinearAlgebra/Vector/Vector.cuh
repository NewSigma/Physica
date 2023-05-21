/*
 * Copyright 2022 WeiBo He.
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
        class Traits<device_obj<Vector<T, Length, MaxLength>>> : public Traits<Vector<T, Length, MaxLength>> {};
    }

    template<class T, size_t Length, size_t MaxLength>
    class device_obj<Vector<T, Length, MaxLength>>
            : public device_obj<ContinuousVector<Vector<T, Length, MaxLength>>>
            , public Utils::device_obj<Utils::Array<T, Length, MaxLength>> {
    public:
        using host_obj = Vector<T, Length, MaxLength>;
    private:
        using Base = device_obj<ContinuousVector<host_obj>>;
        using Storage = Utils::device_obj<Utils::Array<T, Length, MaxLength>>;
    public:
        using Storage::Storage;
        device_obj(const device_obj&);
        device_obj(device_obj&&) noexcept;
        ~device_obj() = default;
        /* Opporators */
        device_obj<Vector<T, Length, MaxLength>> operator=(device_obj<Vector<T, Length, MaxLength>> obj) noexcept;
        using Base::operator=;
        using Storage::operator=;
        using Storage::operator[];
        /* Opporations */
        [[nodiscard]] inline host_obj toHost() const;
        inline void toHost(host_obj& obj);
        using Base::toHost;
        using Storage::resize;
        using Storage::swap;
        /* Getters */
        using Storage::getLength;
        using Storage::data;
    };

    template<class T, size_t Length, size_t MaxLength>
    device_obj<Vector<T, Length, MaxLength>>::device_obj(const device_obj& obj) : Storage(obj) {}

    template<class T, size_t Length, size_t MaxLength>
    device_obj<Vector<T, Length, MaxLength>>::device_obj(device_obj&& obj) noexcept : Storage(std::move(obj)) {}

    template<class T, size_t Length, size_t MaxLength>
    device_obj<Vector<T, Length, MaxLength>>
    device_obj<Vector<T, Length, MaxLength>>::operator=(device_obj<Vector<T, Length, MaxLength>> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class T, size_t Length, size_t MaxLength>
    inline typename device_obj<Vector<T, Length, MaxLength>>::host_obj
    device_obj<Vector<T, Length, MaxLength>>::toHost() const {
        return host_obj(Storage::toHost());
    }

    template<class T, size_t Length, size_t MaxLength>
    inline void device_obj<Vector<T, Length, MaxLength>>::toHost(host_obj& obj) {
        Storage::toHost(obj);
    }

    template<class T, size_t Length, size_t MaxLength, class Allocator>
    inline device_obj<Vector<T, Length, MaxLength, Allocator>> Vector<T, Length, MaxLength, Allocator>::toDevice() const {
        return device_obj<Vector<T, Length, MaxLength>>(*this);
    }

    template<class T, size_t Length, size_t MaxLength, class Allocator>
    inline void Vector<T, Length, MaxLength, Allocator>::toDevice(device_obj<This>& obj) const {
        Storage::toDevice(obj);
    }
}
