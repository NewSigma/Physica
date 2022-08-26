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
#include <unistd.h>
#include <sys/stat.h>
#include <cstring>
#include <cstdio>
#include <cstdarg>
#include <iostream>
#include <wait.h>
#include <cassert>
#include <ctime>
#include "Physica/Utils/Unix/UnixHelper.h"

namespace Physica::Utils {
    void statCheck(const char* __restrict file, struct stat* __restrict buf) {
        if (stat(file, buf) != 0) {
            perror("[Error]: Failed to fetch file stat");
            exit(EXIT_FAILURE);
        }
    }

    int forceRemoveDir(const char* dirPath) {
        DIR* dir = opendir(dirPath);
        if (dir == nullptr)
            return 1;

        dirent* item;
        struct stat st{};
        auto buffer = getPathBuffer();
        while ((item = readdir(dir)) != nullptr) {
            const char* filename = item->d_name;
            if(strcmp(filename, ".") == 0 || strcmp(filename,"..") == 0) {
                continue;
            }

            sprintf(buffer.get(), "%s/%s", dirPath, filename);
            statCheck(buffer.get(), &st);
            if (S_ISDIR(st.st_mode)) {
                if (forceRemoveDir(buffer.get()) != 0)
                    return 1;
            }
            else {
                if(unlink(buffer.get()) != 0)
                    return 1;
            }
        }
        closedir(dir);
        return rmdir(dirPath);
    }
    /**
     * @return a array whose length is enough to store path name.
     */
    std::unique_ptr<char[]> getPathBuffer() {
        const long maxPathLength = pathconf("/", _PC_PATH_MAX);
        return std::unique_ptr<char[]>(new char[maxPathLength]);
    }

    std::unique_ptr<char[]> makePath(const char* format, ...) {
        const long maxPathLength = pathconf("/", _PC_PATH_MAX);
        auto buffer = getPathBuffer();

        std::va_list args;
        va_start(args, format);
        int count = vsnprintf(buffer.get(), maxPathLength, format, args);
        va_end(args);
        if (count == maxPathLength) {
            std::cerr << "[Error]: Working directory is too deep to execute.\n";
            exit(EXIT_FAILURE);
        }
        return buffer;
    }
    /**
     * Ensure the file(regular file or directory) is not exist. If it exists, remove it.
     * @param path
     */
    void ensureNotExist(const char* path) {
        struct stat st{};
        errno = 0;
        if (stat(path, &st) != 0) {
            if (errno != ENOENT) {
                perror("[Error]: Failed to fetch stat of file");
                exit(EXIT_FAILURE);
            }
        }
        if(S_ISDIR(st.st_mode))
            forceRemoveDir(path);
        else
            unlink(path);
    }
    /**
     * Generate runtime folder. If it is not empty, rename it by appending date to dirname.
     */
    void genRuntimeDir() {
        struct stat st{};
        errno = 0;
        if (stat("runtime", &st) != 0) {
            if (errno != ENOENT) {
                perror("[Error]: Failed to fetch stat of runtime.");
                exit(EXIT_FAILURE);
            }
        }
        else if (S_ISDIR(st.st_mode)) {
            time_t now = time(nullptr);
            auto localTime = std::localtime(&now);
            constexpr size_t length = 29;
            char buffer[length];
            [[maybe_unused]] size_t count = strftime(buffer, length, "runtime_%F_%T", localTime); //TODO should be create time rather than rename time
            assert(count != 0);
            errno = 0;
            if (rename("runtime", buffer) != 0) {
                perror("[Error]: Failed to remove old runtime folder.");
                exit(EXIT_FAILURE);
            }
        }
        else
            unlink("runtime");
        mkdir("runtime", S_IRWXU);
    }
    /**
     * Wait for child processes. If anything went wrong, print a error message.
     */
    pid_t waitCheck(const char* message) {
        int status;
        pid_t finished = wait(&status);
        if (finished <= 0) {
            fprintf(stderr, "[Error]: Failed to wait for chile processes.\n");
            exit(EXIT_FAILURE);
        }

        int error = -1;
        if (WIFEXITED(status))
            error = WEXITSTATUS(status);
        if (error != 0) {
            fprintf(stderr, "%s\n", message);
            exit(EXIT_FAILURE);
        }
        return finished;
    }
    /**
     * Copy and sync at once. Because other program will read th results later.
     */
    void copyFile(const char* from, const char* to) {
        FILE* fromFile = fopen(from, "r");
        FILE* toFile = fopen(to, "w");
        assert(fromFile != nullptr && toFile != nullptr && fromFile != toFile);
        int c;
        while ((c = fgetc(fromFile)) != EOF)
            fputc(c, toFile);
        fclose(fromFile);
        fclose(toFile);
    }
    /**
     * Copy all regular files under a dictory to another directory.(Not recursively)
     */
    void copyFileFromDir(const char* fromDirPath, const char* toDirPath) {
        DIR* fromDir = opendir(fromDirPath);
        if (fromDir == nullptr) {
            fprintf(stderr, "[Error]: Failed to open dir: %s", strerror(errno));
            exit(EXIT_FAILURE);
        }

        dirent* item;
        struct stat st{};
        auto fromFileName = getPathBuffer();
        auto toFileName = getPathBuffer();
        while ((item = readdir(fromDir)) != nullptr) {
            const char* filename = item->d_name;
            if(strcmp(filename, ".") == 0 || strcmp(filename,"..") == 0) {
                continue;
            }

            sprintf(fromFileName.get(), "%s/%s", fromDirPath, filename);
            statCheck(fromFileName.get(), &st);
            if (S_ISREG(st.st_mode)) {
                sprintf(toFileName.get(), "%s/%s", toDirPath, filename);
                copyFile(fromFileName.get(), toFileName.get());
            }
        }
        closedir(fromDir);
    }

    void mkdirCheck(const char* path, mode_t mode) {
        int err = mkdir(path, mode);
        if (err != 0) {
            perror("[Error]: Failed to make directory");
            exit(EXIT_FAILURE);
        }
    }
}
