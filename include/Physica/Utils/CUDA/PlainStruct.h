/*
 * Copyright 2023 WeiBo He.
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

namespace Physica {
    template<class Derived> class PlainStruct;
    template<> class PlainStruct<void> {};
    /**
     * \class PlainStruct pass objects to cuda kernel ignoring constructors and destructors because resource control is cpu's duty.
     */
    template<class Derived>
    class alignas(Derived) PlainStruct<const Derived> {
        char anonymous[sizeof(Derived)];
    public:
        [[nodiscard]] __host__ __device__ const Derived& getDerived() const noexcept { return *reinterpret_cast<const Derived*>(this); }
        [[nodiscard]] __host__ __device__ Derived& getConstCastDerived() const noexcept { return *reinterpret_cast<Derived*>(const_cast<PlainStruct*>(this)); }
    };

    template<class Derived>
    class alignas(Derived) PlainStruct : public PlainStruct<const Derived> {
        using Base = PlainStruct<const Derived>;
    public:
        using Base::getDerived;
        [[nodiscard]] __device__ Derived& getDerived() noexcept { return *reinterpret_cast<Derived*>(this); }
    };

    template<class T>
    __host__ __device__ inline PlainStruct<const T> asStruct(const T& obj) noexcept {
        return PlainStruct<const T>(reinterpret_cast<const PlainStruct<const T>&>(obj));
    }

    template<class T>
    __host__ __device__ inline PlainStruct<T> asStruct(T& obj) noexcept {
        return PlainStruct<T>(reinterpret_cast<PlainStruct<T>&>(obj));
    }
}
