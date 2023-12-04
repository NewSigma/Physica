/*
 * Copyright 2023 WeiBo He.
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
#include <iostream>
#include "Physica/Core/MultiPrecision/Differentiable.h"
#include "Physica/Core/MultiPrecision/AutoDiffGuard.h"
#include "Physica/Core/Physics/MD/ForceModel/SilveraGoldman.h"

using namespace Physica::Core;

using PlainScalar = Scalar<Double>;
using ScalarType = Differentiable<PlainScalar, DiffMode::Reverse>;

class PressTest {
    using RandomGenerator = std::mt19937;
    using MDCellType = MDCell<ScalarType>;
    using LatticeMatrix = typename MDCellType::LatticeMatrix;
    using PositionMatrix = typename MDCellType::PositionMatrix;
    using MassVector = typename MDCellType::MassVector;
    constexpr static double mass = PhyConst<AU>::atomMass(1) * 2;
    constexpr static size_t numMolecular = 108;
    constexpr static double pair_cutoff = 15;
public:
    static void run() {
        const AutoDiffGuard<PlainScalar> guard{};
        RandomGenerator gen{};
        const ScalarType volume = 8000;
        const auto cell = makeSystem(gen, volume);
        SilveraGoldman<ScalarType, true> sg(pair_cutoff);
        sg.potentialEnergy(cell).reverse();
        const ScalarType press_diff = -volume.getTangent();
        const ScalarType press = sg.virial(cell).trace() / ScalarType(3);
        if (!scalarNear(press_diff.getValue(), press.getValue(), 1E-15))
            exit(EXIT_FAILURE);
    }
private:
    static MDCellType makeSystem(RandomGenerator& gen, ScalarType volume) {
        LatticeMatrix lattice = MDCellType::LatticeMatrix::unitMatrix(3);
        const ScalarType latticeConst = cbrt(volume);
        lattice *= latticeConst;

        PositionMatrix pos(numMolecular, 3);
        pos.random_uniform(gen);
        pos *= latticeConst;

        MassVector massVec(numMolecular, mass);
        return MDCellType(std::move(lattice), std::move(pos), std::move(massVec));
    }
};

int main() {
    SilveraGoldman<ScalarType, true> sg(1.0);
    {
        const AutoDiffGuard<PlainScalar> guard{};
        ScalarType r = 2.0;
        const ScalarType r2 = square(r);
        sg.pot_functor(0, 0, r, r2).reverse();
        const ScalarType f = -r.getTangent();
        const ScalarType f1 = sg.force_functor(0, 0, r, r2);
        if (!scalarNear(f, f1, 1E-15))
            return 1;
    }
    {
        const AutoDiffGuard<PlainScalar> guard{};
        ScalarType r = 2.0;
        const ScalarType r2 = square(r);
        sg.force_functor(0, 0, r, r2).reverse();
        const ScalarType fc = -r.getTangent();
        const ScalarType fc1 = sg.forceConst_functor(r, r2);
        if (!scalarNear(fc, fc1, 1E-15))
            return 1;
    }
    PressTest::run();
    return 0;
}
