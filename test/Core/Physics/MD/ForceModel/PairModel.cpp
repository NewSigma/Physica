/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.h"
#include "Physica/Core/Physics/MD/ForceModel/LJModel.h"
#include "Physica/Core/Physics/MD/ForceModel/SilveraGoldman.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Test.h"

using namespace Physica;
using dfloat = Diff<float64, DiffMode::Reverse, 1>;

int main() {
    {
        LJModel<dfloat> lj(1.0, 1.0);
        dfloat r = 1.0;
        lj.pot_functor(0, 0, r, square(r)).reverse();
        const float64 f = -r.grad();
        const auto f1 = lj.force_functor(0, 0, r, square(r));
        expect(scalarNear(f.value(), f1.value(), 1E-15));
    }
    return 0;
}
