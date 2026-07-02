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
#include "Physica/CRCoro.h"
#include "Test.h"

using namespace Physica;

namespace {
    struct A : public CRCoro<A> {
        int i;
        A() = default;
        A(int i_) : i(i_) {}
    };

    A fn() {
        const A x{123};
        co_return x; // Regression test for returning a const value
    }
}

int main() {
    expect(fn().i == 123);
}
