/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "AbstractLogger.h"

namespace Physica {
    class PHYSICA_API FileLogger final : public AbstractLogger {
        int fd;
    public:
        FileLogger(const char* filename, bool trunc);
        ~FileLogger() override;
        /* Operations */
        void log(LogBuffer& buffer) override final;
    };
}
