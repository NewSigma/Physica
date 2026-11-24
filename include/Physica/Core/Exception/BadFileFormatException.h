/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "IOException.h"
#include "Physica/Macro.h"

namespace Physica {
    class PHYSICA_API BadFileFormatException : public IOException {
        using This = BadFileFormatException;
    public:
        BadFileFormatException(const char* msg) : IOException(msg) {}
        BadFileFormatException(const This&) = default;
        BadFileFormatException(This&&) noexcept = default;
        ~BadFileFormatException() noexcept override = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
    };
}
