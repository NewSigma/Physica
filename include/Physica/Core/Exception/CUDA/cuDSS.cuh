/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Macro.h"

#ifdef PHYSICA_CUDA
#include <cudss.h>

namespace Physica {
    PHYSICA_API void check(cudssStatus_t err);
}
#endif
