/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <cstdarg>
#include "Physica/Logger/LogBuffer.h"
#include "Physica/Logger/LoggerType.h"

namespace Physica {
    /**
     * Abstract father class for all loggers.
     */
    class PHYSICA_API AbstractLogger {
    public:
        LogLevel localLevel;
    public:
        explicit AbstractLogger(LogLevel level = LogLevel::Global) : localLevel(level) {}
        AbstractLogger(const AbstractLogger& logger) = delete;
        AbstractLogger(AbstractLogger&& logger) noexcept = delete;
        virtual ~AbstractLogger() = default;
        /* Operators */
        AbstractLogger& operator=(const AbstractLogger&) = delete;
        AbstractLogger& operator=(AbstractLogger&&) noexcept = delete;
        /* Operations */
        virtual void log(LogBuffer& buffer) = 0;
        /* Getters */
        [[nodiscard]] LogLevel getCurrentLevel() const noexcept;
    };
}
