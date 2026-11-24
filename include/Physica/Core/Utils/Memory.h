/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <cstddef>
#include "Physica/Macro.h"

namespace Physica {
    PHYSICA_API void memcpy2D(void* dst, size_t dstride, const void* src, size_t sstride, size_t length, size_t numRow) noexcept;
    PHYSICA_API void memset2D(void* dst, size_t dstride, int value, size_t length, size_t numRow) noexcept;
    PHYSICA_API void memcpy_nt(void* dst, const void* src, size_t n) noexcept;
}
