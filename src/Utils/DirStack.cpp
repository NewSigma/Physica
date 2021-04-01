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
#include "Physica/Utils/DirStack.h"
#include <cstring>
#include <cassert>

namespace Physica::Utils {
    DirStack::DirStack(const char* path) {
        char ch = path[0];
        if (ch != '/')
            throw std::invalid_argument("Path must start with '/'\n");
        size_t count = 1;
        size_t dirStart = 1;
        size_t dirEnd = 1;
        do {
            if (ch == '/') {
                if (dirStart == dirEnd)
                    goto ignore;
                cutParh(path, dirStart, dirEnd);
                dirStart = ++dirEnd;
            }
            else {
                ++dirEnd;
            }
        ignore:
            ch = path[count++];
        } while(ch != '\0');

        if (dirStart != dirEnd)
            cutParh(path, dirStart, dirEnd);
    }

    DirStack::~DirStack() {
        for (auto dir : dirs)
            delete[] dir;
    }

    void DirStack::push(const char* dir) {
        const size_t length = strlen(dir);
        char* copy = new char[length + 1];
        memcpy(copy, dir, length + 1);
        dirs.push_back(copy);
    }

    std::unique_ptr<const char[]> DirStack::pop() {
        auto ite = dirs.crbegin();
        const char* value = *ite;
        dirs.pop_back();
        return std::unique_ptr<const char[]>(value);
    }
    
    std::string DirStack::toPath() const {
        std::string result{};
        for (const auto& dir : dirs) {
            result.push_back('/');
            result += dir;
        }
        if (result.empty())
            result.push_back('/');
        return result;
    }

    void DirStack::cutParh(const char* path, size_t startPos, size_t endPos) {
        assert(startPos < endPos);
        const size_t dirLength = endPos - startPos + 1;
        char* __restrict dir = new char[dirLength];
        const size_t dirLength_1 = dirLength - 1;
        memcpy(dir, path + startPos, dirLength_1);
        dir[dirLength_1] = '\0';
        dirs.push_back(dir);
    }
}
