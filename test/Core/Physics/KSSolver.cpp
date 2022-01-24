/*
 * Copyright 2021-2022 WeiBo He.
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
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "Physica/Core/Physics/ElectronicStructure/DFT/KSSolver.h"

using namespace Physica::Core;
using ScalarType = Scalar<Double, false>;

int main() {
    CrystalCell Si({5, 0, 0, 0, 5, 0, 0, 0, 5}, {0.5, 0.5, 0.5}, {14});
    KSSolver<ScalarType> solver = KSSolver<ScalarType>(Si, 0.8, 100, 100, 100);
    solver.solve(1E-3, 100);
    return 0;
}
