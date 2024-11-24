/*
 * Copyright 2021-2024 Weibo He.
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
#include <cassert>
#include <cstring>
#include "Physica/Core/Utils/DirStack.h"

namespace Physica::Core {
    class Test {
    public:
        static std::vector<const char*>& getDirs(DirStack& stack) { return stack.dirs; }
    };
}

int main() {
    using namespace Physica::Core;
    DirStack stack("/home/user/Program");
    auto& dirs = Test::getDirs(stack);
    if (dirs.size() != 3)
        return 1;
    if (strcmp(dirs[0], "home") != 0)
        return 1;
    if (strcmp(dirs[1], "user") != 0)
        return 1;
    if (strcmp(dirs[2], "Program") != 0)
        return 1;

    stack.pop();
    if (stack.toPath() != "/home/user")
        return 1;
    stack.push("AnotherProgram");
    if (stack.toPath() != "/home/user/AnotherProgram")
        return 1;
    return 0;
}
