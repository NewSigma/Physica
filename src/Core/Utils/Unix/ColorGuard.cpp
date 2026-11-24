/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Utils/Unix/ColorGuard.h"

using namespace Physica;

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
    {ALLCOLORS("3", ""), ALLCOLORS("3", "1;")},
    {ALLCOLORS("4", ""), ALLCOLORS("4", "1;")}
};
