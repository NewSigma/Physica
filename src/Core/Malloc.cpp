/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
// Global overloads of operator new/delete, compile it with downstream project if necessary
#ifdef PHYSICA_MIMALLOC
    // Mimalloc includes these headers and they contains always_inline; pre-include them before we attach noinline attributes.
    #include <new>
    #include <mimalloc.h>
    // We are using __builtin_operator_new/delete and noinline is required to keep the builtin semantics throughout the ThinLTO pipeline.
    #pragma clang attribute push(__attribute__((noinline)), apply_to = function)
    #include <mimalloc-new-delete.h>
    #pragma clang attribute pop
#endif
