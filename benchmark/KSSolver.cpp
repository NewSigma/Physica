/*
 * Copyright 2022 WeiBo He.
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
#include "Physica/Utils/TestHelper.h"
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "Physica/Core/Physics/ElectronicStructure/DFT/KSSolver.h"
#include "Physica/Core/Physics/ElectronicStructure/DFT/XCProvider/LDA.h"
#include "Physica/Core/Physics/ElectronicStructure/DFT/BandGrid.h"
#include "Physica/Utils/Cycler.h"

using namespace Physica::Core;
using ScalarType = Scalar<Double, false>;
using ComplexType = ComplexScalar<ScalarType>;

int main() {
    using namespace Physica::Utils;
    Cycler::init();
    CrystalCell Si({5, 0, 0, 0, 5, 0, 0, 0, 5}, {0.5, 0.5, 0.5}, {14});
    ScalarType cutEnergy(0.8);
    {
        BandGrid<ScalarType, false> grid(cutEnergy, Si.reciprocal().getLattice(), 1, 1, 1, 14);
        auto solver = KSSolver<ScalarType, LDA<ScalarType, LDAType::HL, false>>(Si, cutEnergy, 50, std::move(grid));
        const auto from = Cycler::tic();
        solver.solve(1E-3, 100);
        const auto to = Cycler::toc();
        std::cout << "Time use: " << Cycler::toSeconds(to - from) << '\n';
    }
    return 0;
}
