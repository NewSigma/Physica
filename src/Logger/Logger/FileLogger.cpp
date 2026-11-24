/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include <cstring>
#include <iostream>
#include <fcntl.h>
#include <unistd.h>
#include "Physica/Core/Exception/IOException.h"
#include "Physica/Logger/Logger/FileLogger.h"

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
