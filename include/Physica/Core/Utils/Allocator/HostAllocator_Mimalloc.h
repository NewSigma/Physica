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

#include "HostAllocator.h"

namespace Physica {
    template<class T, size_t Align>
    T* HostAllocator<T, Align>::reallocate_mimalloc(T* p, size_t new_size, [[maybe_unused]] size_t old_size) noexcept {
        assert(new_size > 0 && "[Error]: Allocate nothing");
        assert(p != nullptr || old_size == 0); // According to [1], the behavior is well defined now
        void* new_p{};
        if constexpr (OverAlign)
            new_p = mi_realloc_aligned(p, new_size * sizeof(T), Align);
        else
            new_p = realloc(p, new_size * sizeof(T));
        assert(new_p != nullptr);
        return reinterpret_cast<T*>(new_p);
    }
}
