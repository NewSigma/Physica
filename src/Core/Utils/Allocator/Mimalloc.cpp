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
#include <cassert>
#include <mimalloc.h>
#include "Physica/Core/Utils/Allocator/HostAllocator.h"
/**
 * Reference:
 * [1] N3322; https://www.open-std.org/jtc1/sc22/wg14/www/docs/n3322.pdf
 */
void* Physica::Internal::reallocate_mimalloc(void* p, size_t new_size, [[maybe_unused]] size_t old_size, size_t align) noexcept {
    assert(new_size > 0 && "[Error]: Reject bad pattern");
    assert(p != nullptr || old_size == 0); // According to [1], the behavior is well defined now
    void* new_p{};
    if (align == 0)
        new_p = mi_realloc(p, new_size);
    else
        new_p = mi_realloc_aligned(p, new_size, align);
    assert(new_p != nullptr);
    return new_p;
}
