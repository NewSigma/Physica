/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Scalar/Diff.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using dfloat = Diff<float64, DiffMode::Forward, 2>;

namespace {
    void tupleTest() {
        auto [p1, p2, p3] = ScalarPtr<dfloat>{};
        static_assert(std::same_as<decltype(p1), T*>);
        static_assert(std::same_as<decltype(p3), T*>);
    }
}

int main() {
    static_assert(std::random_access_iterator<ScalarPtr<Diff<float64, DiffMode::Forward>>>);
    static_assert(std::formattable<ScalarRef<Diff<float64, DiffMode::Forward>>, char>);
    tupleTest();
    return 0;
}
