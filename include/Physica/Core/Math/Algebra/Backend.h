/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

namespace Physica {
    enum class Backend : char {
        Dispatch = 0,
        Base,
        MKL,
        cuBLAS,
    };
}
