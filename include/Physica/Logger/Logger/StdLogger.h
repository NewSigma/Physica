/*
 * Copyright 2020-2021 WeiBo He.
 *
 * This file is part of Physica.
 *
 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include "AbstractLogger.h"

namespace Physica::Logger {
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
