/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2024 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Scalar/Integer.h"

namespace Physica {
    Integer factorial(const Integer& i) {
        Integer result = 1;
        Integer temp = 1;
        while(temp < i)
            result *= ++temp;
        return result;
    }
}
