/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2024 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <stdexcept>
#include "Physica/Macro.h"

namespace Physica {
    class PHYSICA_API BadConvergenceException : public std::runtime_error {
    public:
        BadConvergenceException(const char* msg_) noexcept : std::runtime_error(msg_) {}
    };
}
