/*
 * Copyright 2024 Weibo He.
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

#include <ostream>
#include "Physica/Macro.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] LLVM; https://github.com/llvm/llvm-project/blob/main/clang/include/clang/AST/ASTDumperUtils.h
     */
    class PHYSICA_API ColorGuard {
    public:
        enum class Color {
            Black = 0,
            Red,
            Green,
            Yellow,
            Blue,
            Magenta,
            Cyan,
            White,
            SavedColor,
            Reset
        };

        static const char colorcodes[2][2][8][10];
    private:
        std::ostream& os;
    public:
        inline ColorGuard(std::ostream& os_, Color color, bool bold, bool background = false);
        inline ~ColorGuard();
    };

    inline ColorGuard::ColorGuard(std::ostream& os_, Color color, bool bold, bool background) : os(os_) {
        os << colorcodes[background ? 1 : 0][bold ? 1 : 0][int(color) & 7];
    }

    inline ColorGuard::~ColorGuard() {
        os << "\033[0m";
    }
}
