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
#include <fstream>
#include <gperftools/profiler.h>
#include "Physica/Core/IO/Poscar.h"
#include "Physica/Core/Physics/MD/EnergyMinimizer.h"
#include "Physica/Core/Physics/MD/ForceModel/Q_TIP4P.h"
#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"

using namespace Physica::Core;
using ScalarType = Scalar<Double, false>;
using Optimizer = SteepestDescent<ScalarType, Dynamic>;
using Minimizer = EnergyMinimizer<ScalarType, ScalarType>;
using MDCellType = typename Minimizer::MDCellType;
using ForceModel = Q_TIP4P<ScalarType, ScalarType>;

int main() {
    std::ifstream fin("POSCAR");
    Poscar poscar{};
    fin >> poscar;

    MDCellType cell(poscar);
    cell.scale(PhyConst<AU>::angstormToBohr(1));
    auto model = ForceModel(cell, PhyConst<AU>::angstormToBohr(9));
    ForceModel::sortPosition(cell);
    Optimizer sd(1E-5);
    auto minimizer = Minimizer(cell);
    minimizer.init(model, sd);
    size_t step = 0;
    while (!minimizer.pos_step<ForceModel, SequentialExecutor, Optimizer>(model, sd)) {
        step += 1;
        if (step == 100)
            break;
    }

    std::ofstream fout("POSCAR");
    auto optimized_cell = minimizer.getCell();
    optimized_cell.scale(PhyConst<AU>::bohrToAngstorm(1));
    fout << Poscar(optimized_cell, poscar.getElementTypes(), poscar.getNumOfEachType());
    return 0;
}
