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
#include <fstream>
#include "Physica/Core/IO/Poscar.h"
#include "Physica/Core/Physics/MD/EnergyMinimizer.h"
#include "Physica/Core/Physics/MD/ForceModel/Q_TIP4P.h"
#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"
#include "Physica/Utils/Unix/TempFile.h"

using namespace Physica::Core;
using ScalarType = Scalar<Double>;
using Optimizer = SteepestDescent<ScalarType, Dynamic>;
using Minimizer = EnergyMinimizer<ScalarType, ScalarType>;
using MDCellType = typename Minimizer::MDCellType;
using ForceModel = Q_TIP4P<ScalarType, ScalarType>;

const static char* data1 = "\n"
                           "1.0\n"
                           "  -6.362084866    6.296391964    17.81582451\n"
                           " -0.1034414023    6.610432148 -0.04331437126\n"
                           "  -2.796585083 -0.04086640105 -0.01148547977\n"
                           " H O\n"
                           " 8 4\n"
                           "C\n"
                           "-1.161022305  3.972291231  2.481887341\n"
                           "-2.088842869  2.777185917   2.64332962\n"
                           "-1.152102232  2.639525414 0.7057020664\n"
                           "-8.427772522  8.855648994  17.30949211\n"
                           "-2.992442846  5.872242451  2.691781044\n"
                           "-1.11223042  5.951196671  1.491070032\n"
                           "-7.404219627  6.947033882  17.55684662\n"
                           "-6.544594765  12.38838673   17.3673954\n"
                           "-1.181259394  3.005724192  2.426948071\n"
                           "-7.517476559   8.71113205  17.57943535\n"
                           "-1.101659656  5.738564014   2.43592453\n"
                           "-7.258571434  12.59273052  17.58112335\n";

Poscar makePoscar() {
    auto tmp = Physica::Utils::TempFile("/tmp/tmpXXXXXX");
    {
        std::ofstream fout(tmp.getName());
        fout << data1;
    }
    std::ifstream fin(tmp.getName());
    Poscar poscar{};
    fin >> poscar;
    return poscar;
}

int main() {
    auto poscar = makePoscar();

    MDCellType cell(poscar);
    cell.scale(PhyConst<AU>::angstormToBohr(1));
    auto model = ForceModel(cell, PhyConst<AU>::angstormToBohr(9));
    ForceModel::sortPosition(cell);
    Optimizer sd(1, 1E-4, 0.1);
    auto minimizer = Minimizer(cell);
    minimizer.pre_pos_step<ForceModel, SequentialExecutor, Optimizer>(model, sd);
    for (int i = 0; i < 100; ++i) {
        minimizer.pos_step<ForceModel, SequentialExecutor, Optimizer>(model, sd);
    }

    auto optimized_cell = minimizer.getCell();
    optimized_cell.scale(PhyConst<AU>::bohrToAngstorm(1));
    std::cout << Poscar(optimized_cell, poscar.getElementTypes(), poscar.getNumOfEachType()) << '\n';
    return 0;
}
