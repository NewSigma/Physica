/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Macro.h"
#ifdef PHYSICA_MKL
    #include <mkl_vsl.h>
#endif

namespace Physica {
#ifdef PHYSICA_MKL
    PHYSICA_API void check_vsl(int err);
#else
    #define check_vsl(expr)
#endif
}
