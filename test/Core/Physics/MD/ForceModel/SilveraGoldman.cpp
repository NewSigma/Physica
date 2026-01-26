/*
 * Copyright 2023-2026 Weibo He.
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
#include "Physica/Core/Physics/MD/ForceModel/SilveraGoldman.h"
#include "Test.h"

using namespace Physica;

using Tv = float64;
using dfloat = Diff<Tv, DiffMode::Reverse, 1>;

// FIXME: Re-enable this test
constexpr bool Disable = sizeof(int) == 0;

class PressTest {
    using RandomSource = Random<>;
    using MDCellType = MDCell<dfloat>;
    using LatticeMatrix = MDCellType::LatticeMatrix;
    using PositionMatrix = MDCellType::PositionMatrix;
    using MassVector = MDCellType::MassVector;
    constexpr static double mass = PhyConst<AU>::atomMass(1) * 2;
    constexpr static size_t numMolecular = 108;
    constexpr static double pair_cutoff = 15;
public:
    static void run() {
        if constexpr (Disable) {
            dfloat volume = 8000;
            const auto cell = makeSystem(volume);
            SilveraGoldman<dfloat, true> sg(pair_cutoff);
            sg.potentialV(cell).reverse();
            const float64 press_diff = -volume.grad();
            const float64 press = sg.virial(cell).trace().value() / float64(3);
            expect(scalarNear(press_diff, press, 1E-14));
        }
    }
private:
    static MDCellType makeSystem(dfloat& volume) {
        if constexpr (Disable) {
            LatticeMatrix lattice = MDCellType::LatticeMatrix::identity(3);
            const auto latticeConst = cbrt(volume);
            lattice *= latticeConst.value();

            PositionMatrix pos(numMolecular, 3);
            pos.random_uniform<RandomSource>();
            pos *= latticeConst.value();

            MassVector massVec(numMolecular, mass);
            return MDCellType(std::move(lattice), std::move(pos), std::move(massVec));
        }
        return {};
    }
};

int main() {
    if constexpr (Disable) {
        SilveraGoldman<dfloat, true> sg(1.0);
        {
            dfloat r = 2.0;
            const auto r2 = square(r);
            sg.pot_functor(0, 0, r, r2).reverse();
            const Tv f = -r.grad();
            const Tv f1 = sg.force_functor(0, 0, r, r2).value();
            expect(scalarNear(f, f1, 1E-15));
        }
        {
            dfloat r = 2.0;
            const auto r2 = square(r);
            sg.force_functor(0, 0, r, r2).reverse();
            const Tv fc = -r.grad();
            const Tv fc1 = sg.forceConst_functor(r, r2).value();
            expect(scalarNear(fc, fc1, 1E-15));
        }
    }
    PressTest::run();
    return 0;
}
