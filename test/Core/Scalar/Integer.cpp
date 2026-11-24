/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Scalar/Integer.h"
#include "Test.h"

using namespace Physica;

static_assert(factorial(0) == 1);
static_assert(factorial(20) == 2432902008176640000ULL);
static_assert(doubleFactorial(16) == 10321920ULL);
static_assert(doubleFactorial(17) == 34459425ULL);
static_assert(doubleFactorial(33) == 6332659870762850625ULL);

int main() {
    expect(factorial(Integer(5)) == Integer(120));
    return 0;
}
