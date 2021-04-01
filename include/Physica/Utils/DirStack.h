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
#include <vector>
#include <memory>

namespace Physica::Utils {
    class DirStack {
        std::vector<const char*> dirs;
    public:
        DirStack(const char* path);
        ~DirStack();
        /* Operations */
        void push(const char* dir);
        std::unique_ptr<const char[]> pop();
        std::string toPath() const;
    private:
        void cutParh(const char* path, size_t startPos, size_t endPos);
        friend class Test;
    };
}
