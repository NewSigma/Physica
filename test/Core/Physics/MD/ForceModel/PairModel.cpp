/*
 * Copyright 2023-2025 Weibo He.
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
