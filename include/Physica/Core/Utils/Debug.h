/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Macro.h"

namespace Physica {
    template <class... T>
    __host__ __device__ consteval void print_type() {
        static_assert(sizeof...(T) < 0, "[Info]: Printing type T");
    }

    template <class... T>
    __host__ __device__ consteval void print_type(T&&...) {
        static_assert(sizeof...(T) < 0, "[Info]: Printing type T");
    }
}
