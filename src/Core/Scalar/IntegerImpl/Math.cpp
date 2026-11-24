/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2024 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
