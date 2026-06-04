/*
 * Copyright 2023-2026 Weibo He.
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
#include "Macro.h"

namespace Physica {
    /**
     * \class PlainStruct pass objects to cuda kernel ignoring constructors and destructors because resource control is cpu's duty.
     */
    template<class T>
    class alignas(T) PlainStruct : public PlainStruct<const T> {
        using Base = PlainStruct<const T>;
    public:
        using Base::getDerived;
        [[nodiscard, gnu::always_inline, gnu::nodebug]] __host__ __device__ constexpr T& getDerived() noexcept { return *reinterpret_cast<T*>(this); }
    };

    template<class T>
    class alignas(T) PlainStruct<const T> {
        using This = PlainStruct<const T>;
    public:
        using Derived = T;
    private:
        std::array<std::byte, sizeof(T)> bytes;
    public:
        [[nodiscard, gnu::always_inline, gnu::nodebug]] __host__ __device__ constexpr const T& getDerived() const noexcept { return *reinterpret_cast<const T*>(this); }
        [[nodiscard, gnu::always_inline, gnu::nodebug]] __host__ __device__ constexpr T& getConstCastDerived() const noexcept { return *reinterpret_cast<T*>(const_cast<PlainStruct*>(this)); }
    };

    template<class T>
    [[nodiscard, gnu::always_inline, gnu::nodebug]] __host__ __device__ auto asStruct(const T& obj) noexcept {
        return PlainStruct<const T>(reinterpret_cast<const PlainStruct<const T>&>(obj));
    }

    template<class T>
    [[nodiscard, gnu::always_inline, gnu::nodebug]] __host__ __device__ auto asStruct(T& obj) noexcept {
        return PlainStruct<T>(reinterpret_cast<PlainStruct<T>&>(obj));
    }
}
