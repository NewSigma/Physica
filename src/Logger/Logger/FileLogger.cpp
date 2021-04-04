/*
 * Copyright 2021 WeiBo He.
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
#include <fcntl.h>
#include <unistd.h>
#include <cstring>
#include "Physica/Core/Exception/IOException.h"
#include "Physica/Logger/Logger/FileLogger.h"

namespace Physica::Logger {
    FileLogger::FileLogger(const char* filename)
            : fd(open(filename, O_WRONLY | O_TRUNC | O_CREAT, S_IRUSR | S_IWUSR)) {
        if (fd == -1)
            throw Core::IOException();
    }

    FileLogger::~FileLogger() {
        close(fd);
    }

    void FileLogger::log(LogBuffer& buffer) {
        const std::string msg = buffer.makeMsgString();
        const char* c_str = msg.c_str();
        write(fd, c_str, strlen(c_str) + 1);
        write(fd, "\n", 1);
    }
}