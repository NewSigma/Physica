/*
 * Copyright 2020 WeiBo He.
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
#ifndef PHYSICA_LOGGERRUNTIMEIMPL_H
#define PHYSICA_LOGGERRUNTIMEIMPL_H

#include "Physica/Utils/Cycler.h"

namespace Physica::Logger {
    template<typename T1, typename... Ts>
    inline void writeArgs(T1 head, Ts... args) {
        LoggerRuntime::getInstance().getBuffer().write(head);
        writeArgs(args...);
    }

    inline void writeArgs() { /* Do nothing */ }

    template<typename... Ts>
    void log(size_t logID, Ts... args) {
        size_t time = Utils::Cycler::now();
        writeArgs(logID);
        writeArgs(time);
        writeArgs(args...);
    }
}

#endif
