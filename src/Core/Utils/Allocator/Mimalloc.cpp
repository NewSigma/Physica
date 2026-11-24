/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
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
