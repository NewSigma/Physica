/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Macro.h"

#ifdef PHYSICA_CUDA
#include <cublas_v2.h>

namespace Physica {
    PHYSICA_API void check(cublasStatus_t err);
}
#endif
