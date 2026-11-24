/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <ostream>
#include "Physica/Macro.h"

namespace Physica {
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
