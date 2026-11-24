/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <limits>
#include "Physica/Core/Scalar/Real.h"
#include "Test.h"

using namespace Physica;
using T = float16;

namespace {
    void subnormal() {
        expect(T(0).isSubNormal());
        expect(T(cuda::std::numeric_limits<float16>::denorm_min()).isSubNormal());
        expect(!T(1).isSubNormal());
        expect(!T(-1).isSubNormal());
        expect(!T(std::numeric_limits<float>::infinity()).isSubNormal());
        expect(!T(-std::numeric_limits<float>::infinity()).isSubNormal());
        expect(!T(std::numeric_limits<float>::quiet_NaN()).isSubNormal());
    }

    void repinf() {
        const T inf = std::numeric_limits<float>::infinity();
        expect(reciprocal(inf).isZero());
        expect(reciprocal(-inf).isZero());
    }
}

int main() {
    subnormal();
    repinf();
    return 0;
}
