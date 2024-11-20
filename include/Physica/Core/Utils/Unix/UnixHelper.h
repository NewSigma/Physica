/*
 * Copyright 2022-2024 Weibo He.
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

#ifdef __linux__

#include <memory>
#include <dirent.h>
#include <sys/stat.h>
#include "Physica/Core/Utils/Container/Array.h"

namespace Physica::Core {
    PHYSICA_API void statCheck(const char* file, struct ::stat* buf);
    PHYSICA_API int forceRemoveDir(const char* dirPath);
    PHYSICA_API long getMaxPathLength();
    PHYSICA_API Array<char> getPathBuffer();
    PHYSICA_API Array<char> makePath(const char* format, ...);
    PHYSICA_API bool fileExists(const char* path);
    PHYSICA_API void ensureNotExist(const char* path);
    PHYSICA_API pid_t waitCheck(const char* message);
    PHYSICA_API void copyFile(const char* from, const char* to);
    PHYSICA_API void copyFileFromDir(const char* fromDirPath, const char* toDirPath);
    PHYSICA_API void mkdirCheck(const char* path, mode_t mode);
}

#endif
