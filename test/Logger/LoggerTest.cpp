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
#include "Physica/PhysicaInit.h"
#include "Physica/Logger/AbstractLogger.h"
#include "Physica/Logger/LoggerRuntime.h"

using namespace Physica::Logger;

int main(int argc, char** argv) {
    (void)argc;
    (void)argv;

    initPhysica();
    auto& logger = getStdoutLogger();

    Info(logger, "Test begin.");
    Warning(logger, "This is %c %s%c", 'a', "Logger", '.');

    char str[] = "This is a dynamic string.";
    Info(logger, "%s", str);

    logger.localLevel = LogLevel::Debug;
    Debug(logger, "This is debug mode.");
    logger.localLevel = LogLevel::Info;
    Debug(logger, "This message should not appear.");
    return 0;
}