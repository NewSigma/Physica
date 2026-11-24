/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2021 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <iosfwd>
#include "AbstractLogger.h"

namespace Physica {
    /**
     * The standard logger that links to stdout and stderr,
     * will be created when LoggerRuntime is initialized.
     */
    class StdLogger final : public AbstractLogger {
        std::ostream& os;
    public:
        StdLogger() = delete;
        ~StdLogger() override = default;
        /* Operations */
        void log(LogBuffer& buffer) override final;
    protected:
        //StdLogger can be created by LoggerRuntime only.
        explicit StdLogger(std::ostream& stream);

        friend class LoggerRuntime;
    };
}
