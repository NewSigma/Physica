/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include <cerrno>
#include <system_error>
#include "Physica/Core/Exception/Syscall.h"
#include "SystemError.h"

using namespace Physica;

void Physica::check_syscall(bool failed) {
    if (failed) [[unlikely]]
        throw SystemError(errno, std::system_category());
}
