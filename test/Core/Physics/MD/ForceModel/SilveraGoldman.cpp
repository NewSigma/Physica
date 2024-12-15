/*
 * Copyright 2023 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.h"
#include "Physica/Core/Physics/MD/ForceModel/SilveraGoldman.h"

using namespace Physica::Core;

using ValueType = float64;
using ScalarType = Diff<ValueType, DiffMode::Reverse, 1>;

class PressTest {
    using RandomType = Random<MT19937>;
    using MDCellType = MDCell<ScalarType>;
    using LatticeMatrix = MDCellType::LatticeMatrix;
    using PositionMatrix = MDCellType::PositionMatrix;
    using MassVector = MDCellType::MassVector;
    constexpr static double mass = PhyConst<AU>::atomMass(1) * 2;
    constexpr static size_t numMolecular = 108;
    constexpr static double pair_cutoff = 15;
public:
    static void run() {
        const ScalarType volume = 8000;
        const auto cell = makeSystem(volume);
        SilveraGoldman<ScalarType, true> sg(pair_cutoff);
        sg.potentialV(cell).reverse();
        const ScalarType press_diff = -volume.grad();
        const ScalarType press = sg.virial(cell).trace() / ScalarType(3);
        if (!scalarNear(press_diff.value(), press.value(), 1E-14))
            exit(EXIT_FAILURE);
    }
private:
    static MDCellType makeSystem(ScalarType volume) {
        LatticeMatrix lattice = MDCellType::LatticeMatrix::unitMatrix(3);
        const ScalarType latticeConst = cbrt(volume);
        lattice *= latticeConst;

        PositionMatrix pos(numMolecular, 3);
        pos.random_uniform<RandomType>();
        pos *= latticeConst;

        MassVector massVec(numMolecular, mass);
        return MDCellType(std::move(lattice), std::move(pos), std::move(massVec));
    }
};

int main() {
    SilveraGoldman<ScalarType, true> sg(1.0);
    {
        ScalarType r = 2.0;
        const ScalarType r2 = square(r);
        sg.pot_functor(0, 0, r, r2).reverse();
        const ScalarType f = -r.grad();
        const ScalarType f1 = sg.force_functor(0, 0, r, r2);
        if (!scalarNear(f, f1, 1E-15))
            return 1;
    }
    {
        ScalarType r = 2.0;
        const ScalarType r2 = square(r);
        sg.force_functor(0, 0, r, r2).reverse();
        const ScalarType fc = -r.grad();
        const ScalarType fc1 = sg.forceConst_functor(r, r2);
        if (!scalarNear(fc, fc1, 1E-15))
            return 1;
    }
    PressTest::run();
    return 0;
}
