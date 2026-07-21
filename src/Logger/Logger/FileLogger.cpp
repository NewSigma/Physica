/*
 * Copyright 2021-2026 Weibo He.
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
module;

#include <cstring>
#include <iostream>
#include <fcntl.h>
#include <unistd.h>
#include "Physica/Core/Exception/IOException.h"

module Physica.Logger.FileLogger;

using namespace Physica;

FileLogger::FileLogger(const char* filename, bool trunc)
        : fd(open(filename, trunc ? (O_WRONLY | O_TRUNC | O_CREAT) : (O_WRONLY | O_APPEND | O_CREAT), S_IRUSR | S_IWUSR)) {
    if (fd == -1)
        throw IOException("[Error]: Failed to open file");
}

FileLogger::~FileLogger() {
    close(fd);
}

void FileLogger::log(LogBuffer& buffer) {
    const std::string msg = buffer.makeMsgString();
    const char* c_str = msg.c_str();
    if (write(fd, c_str, strlen(c_str)) == -1) {
        std::cerr << "[Error]: Failed to write to log file\n";
        exit(EXIT_FAILURE);
    }
    if (write(fd, "\n", 1) == -1) {
        std::cerr << "[Error]: Failed to write to log file\n";
        exit(EXIT_FAILURE);
    }
}
