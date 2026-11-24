/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Logger/LoggerRuntime.h"

namespace Physica {
    LogLevel AbstractLogger::getCurrentLevel() const noexcept {
        return localLevel == LogLevel::Global ? LoggerRuntime::globalLevel : localLevel;
    }
}
