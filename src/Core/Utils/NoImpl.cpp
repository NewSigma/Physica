/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Utils/NoImpl.h"
#include <cassert>
#include <print>

namespace Physica {
    __host__ __device__ void noImpl(const char* msg) noexcept {
        assert(msg != nullptr);
    #ifdef __CUDA_ARCH__
        printf("[Error]: Not implemented\n%s\n", msg);
        __trap();
    #else
        std::println("[Error]: Not implemented\n{}", msg);
        std::abort();
    #endif
    }

    __host__ __device__ void noImpl(std::source_location loc) noexcept {
        noImpl(std::format("File: {}:{}:{}\nFunc: {}", loc.file_name(), loc.line(), loc.column(), loc.function_name()).c_str());
    }
}
