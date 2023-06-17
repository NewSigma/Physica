/*
 * Copyright 2022 WeiBo He.
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

#include <memory>
#include <dirent.h>
#include "Physica/Utils/Container/Array/Array.h"

namespace Physica::Utils {
    void statCheck(const char* __restrict file, struct stat* __restrict buf);
    int forceRemoveDir(const char* dirPath);
    long getMaxPathLength();
    Array<char> getPathBuffer();
    Array<char> makePath(const char* format, ...);
    bool fileExists(const char* path);
    void ensureNotExist(const char* path);
    pid_t waitCheck(const char* message);
    void copyFile(const char* from, const char* to);
    void copyFileFromDir(const char* fromDirPath, const char* toDirPath);
    void mkdirCheck(const char* path, mode_t mode);
}
