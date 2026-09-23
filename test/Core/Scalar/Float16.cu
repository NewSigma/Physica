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
