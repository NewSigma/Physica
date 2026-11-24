/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <stdexcept>
#include <string>
#include "Physica/Macro.h"

namespace Physica {
    class PHYSICA_API IOException : public std::runtime_error {
    public:
        IOException(const char* msg) : std::runtime_error(msg) {}
        IOException(const std::string& msg) : std::runtime_error(msg) {}
    };
}
