/*
 * Copyright 2023-2024 Weibo He.
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
#include "Physica/Core/Physics/MD/ForceModel/LJModel.h"
#include "Physica/Core/Physics/MD/ForceModel/SilveraGoldman.h"
#include "Physica/Core/Math/Random/Random.h"

using namespace Physica;
using dfloat = Diff<float64, DiffMode::Reverse, 1>;
/**
 * Params referenced from [1]
 * 
 * Reference:
 * [1] J. Chem. Phys. 122, 184503 (2005); https://doi.org/10.1063/1.1893956
 */
class ForceConstTest {
    using MDCellType = MDCell<dfloat>;
    using RandomSource = Random<MT19937, 12345>;
    constexpr static unsigned int numMolecular = 32;
    constexpr static double pair_cutoff = 15;
    constexpr static double molarVolume = 31.7;
    constexpr static double mass = PhyConst<AU>::atomMass(1) * 2;
public:
    static void run() {
        SilveraGoldman<dfloat, true> sg(pair_cutoff);
        const auto cell = makeSystem();
        const auto fc = sg.forceConst(cell);
        for (size_t i = 0; i < cell.getDOF(); ++i) {
            for (size_t j = 0; j < cell.getDOF(); ++j) {
                const auto fc1 = sg.forceConst(cell, i, j);
                if (!scalarNear<dfloat>(fc(i, j), fc1, 1E-15))
                    exit(EXIT_FAILURE);
            }
        }
    }
private:
    static MDCellType makeSystem() {
        auto lattice = MDCellType::LatticeMatrix::unitMatrix(3);
        auto pos = MDCellType::PositionMatrix::random_uniform<RandomSource>(numMolecular, 3);
        MDCellType::MassVector massVec(numMolecular, mass);
        MDCellType cell(std::move(lattice), std::move(pos), std::move(massVec));

        const double factor = (std::cbrt(numMolecular * molarVolume / PhyConst<SI>::avogadroNa) / 100) / PhyConst<SI>::bohrRadius;
        cell.scale(factor);
        return cell;
    }
};

int main() {
    {
        LJModel<dfloat> lj(1.0, 1.0);
        dfloat r = 1.0;
        lj.pot_functor(0, 0, r, square(r)).reverse();
        const float64 f = -r.grad();
        const auto f1 = lj.force_functor(0, 0, r, square(r));
        if (!scalarNear(f.value(), f1.value(), 1E-15))
            return 1;
    }
    ForceConstTest::run();
    return 0;
}
