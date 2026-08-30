/*
 * Copyright 2026 Weibo He.
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
#include <algorithm>
#include <bit>
#include <cassert>
#include <cstdint>
#include "Physica/Core/Utils/MetaProgramming.h"

namespace Physica {
    enum class HandleType : int8_t {
        MPI_Comm,
        MPI_Dtype,
    };
    /**
     * \class Handle: A zero-cost wrapper for opaque handles in 3rdparty C libraries, dedicated to provide:
     * 1. Type safety
     * 2. ABI stability: for example, OpenMPI implement handles using void*, while MPICH uses int
     */
    template<HandleType HT>
    class Handle {
        using This = Handle;
        constexpr static size_t SizeLimit = sizeof(void*);

        std::array<std::byte, SizeLimit> opaque{};
    public:
        constexpr Handle() = default;
        constexpr explicit Handle(auto handle) noexcept;
        constexpr Handle(const This&) = default;
        constexpr Handle(This&&) noexcept = default;
        constexpr ~Handle() = default;
        /* Operators */
        constexpr This& operator=(This obj) noexcept;
        [[nodiscard]] constexpr bool operator==(const This& other) const noexcept;
        [[nodiscard]] constexpr bool operator!=(const This& other) const noexcept = default;

        template<class T>
        [[nodiscard]] constexpr explicit operator T() const noexcept;
        /* Operations */
        constexpr void swap(This& obj) noexcept;
    private:
        template<class T>
        consteval static void check() noexcept;
    };

    template<HandleType HT>
    constexpr Handle<HT>::Handle(auto handle) noexcept {
        check<decltype(handle)>();
        opaque.fill(std::byte(0));
        std::ranges::copy(std::bit_cast<std::array<std::byte, sizeof(handle)>>(handle), opaque.begin());
    }

    template<HandleType HT>
    constexpr auto Handle<HT>::operator=(This obj) noexcept -> This& {
        swap(obj);
        return *this;
    }

    template<HandleType HT>
    constexpr bool Handle<HT>::operator==(const This& other) const noexcept {
        return opaque == other.opaque;
    }

    template<HandleType HT>
    template<class T>
    constexpr Handle<HT>::operator T() const noexcept {
        check<T>();
        std::array<std::byte, sizeof(T)> buffer;
        std::ranges::copy(opaque.begin(), opaque.begin() + sizeof(T), buffer.begin());
        return std::bit_cast<T>(buffer);
    }

    template<HandleType HT>
    constexpr void Handle<HT>::swap(This& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        opaque.swap(obj.opaque);
    }

    template<HandleType HT>
    template<class T>
    consteval void Handle<HT>::check() noexcept {
        static_assert(sizeof(T) <= SizeLimit && std::is_trivially_copyable_v<T>);
        static_assert(!instanceof_x<T, Handle>);
    }
}
