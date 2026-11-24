/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Logger/Logger/AbstractLogger.h"
#include "Physica/Logger/LoggerRuntime.h"
#include <unistd.h>

using namespace Physica;

int main() {
    Info(STDOUT_FILENO, "Test begin.");
    Warn(STDERR_FILENO, "This is %c %s%c", 'a', "Logger", '.');

    char str[] = "This is a dynamic string.";
    Info(STDOUT_FILENO, "%s", str);

    AbstractLogger& stdoutLogger = getLogger(1);
    stdoutLogger.localLevel = LogLevel::Debug;
    Debug(STDOUT_FILENO, "This is debug mode.");
    stdoutLogger.localLevel = LogLevel::Info;
    Debug(STDERR_FILENO, "This message should not appear.");
    return 0;
}
