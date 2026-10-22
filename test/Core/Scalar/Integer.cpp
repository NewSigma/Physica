/*
 * Copyright 2026 Weibo He.
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
