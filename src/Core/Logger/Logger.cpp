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
#include "Physica/Logger/Logger.h"

namespace Physica::Logger {
    void Logger::log(const LogInfo* info, ...) {
        printf("[time] [%s]", info->level);
        std::va_list list;
        va_start(list, reinterpret_cast<const void*>(info));
        vprintf(info->format, list);
        va_end(list);
        printf("\n");
    }
}