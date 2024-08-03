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
#include "Physica/Utils/Unix/ColorGuard.h"

namespace Physica::Utils {
    #define COLOR(FGBG, CODE, BOLD) "\033[0;" BOLD FGBG CODE "m"
    #define ALLCOLORS(FGBG, BOLD) {\
        COLOR(FGBG, "0", BOLD),\
        COLOR(FGBG, "1", BOLD),\
        COLOR(FGBG, "2", BOLD),\
        COLOR(FGBG, "3", BOLD),\
        COLOR(FGBG, "4", BOLD),\
        COLOR(FGBG, "5", BOLD),\
        COLOR(FGBG, "6", BOLD),\
        COLOR(FGBG, "7", BOLD)\
    }

    const char ColorGuard::colorcodes[2][2][8][10] = {
        { ALLCOLORS("3", ""), ALLCOLORS("3", "1;") },
        { ALLCOLORS("4", ""), ALLCOLORS("4", "1;") }
    };
}
