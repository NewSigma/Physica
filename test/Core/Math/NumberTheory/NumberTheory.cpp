/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/NumberTheory/NumberTheory.h"
#include "Test.h"

using namespace Physica;

int main() {
    {
        const int i = lcm<int, false>(6, 2);
        expect(i == 6);
    }
    {
        const int i = lcm<int, false>(6, 6);
        expect(i == 6);
    }
    return 0;
}
