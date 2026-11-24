/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Macro.h"
#include <source_location>

namespace Physica {
    [[noreturn, gnu::cold, gnu::nodebug]] PHYSICA_API __host__ __device__ void noImpl(const char* msg) noexcept;
    [[noreturn, gnu::cold, gnu::nodebug]] PHYSICA_API __host__ __device__ void noImpl(std::source_location loc = std::source_location::current()) noexcept;
}
