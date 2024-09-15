/*
 * Copyright 2021-2024 Weibo He.
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
#include "Physica/Core/Math/Random/RandomSeed.h"
#include "Physica/Core/Physics/SolidState/CrystalCell.h"
#include "Physica/Core/Physics/ElectronicStructure/DFT/KSSolver.h"
#include "Physica/Core/Physics/ElectronicStructure/DFT/XCProvider/LDA.h"
#include "Physica/Core/Physics/ElectronicStructure/DFT/BandGrid.h"
#include <Physica/Core/Utils/Cycler.h>

using namespace Physica::Core;
using ScalarType = float64;
using ComplexType = ComplexScalar<ScalarType>;
constexpr bool isSpinPolarized = false;

int main() {
    CrystalCell<ScalarType> Si({{5, 0, 0, 0, 5, 0, 0, 0, 5}, {0.5, 0.5, 0.5}, CrystalCell<ScalarType>::Type::Direct}, {2});
    ScalarType cutEnergy(4.9);

    std::mt19937::result_type seed;
    RandomSeed::rdrand(seed);
    std::mt19937 gen{};
    {
        BandGrid<ScalarType, isSpinPolarized> grid(Si.makeRepLattice(), 1, 1, 1, 2, 4);
        using SolverType = KSSolver<ScalarType, LDA<ScalarType, LDAType::HL, isSpinPolarized>>;
        SolverType solver = SolverType(Si, cutEnergy, 2, std::move(grid), 2);

        {
            using DensityType = typename SolverType::DensityType;
            DensityType last_rho{};
            std::ifstream fin("rho");
            fin >> last_rho;
            solver.getDensity().initDensity(last_rho);
        }
        solver.getOrbits()[0][SpinState::Up].random_normal(gen);

        const auto from = Cycler::tic();
        solver.solve(1E-8, 100);
        const auto to = Cycler::toc();
        std::cout << "Time use: " << Cycler::toSeconds(to - from) << '\n';
        const auto& band = solver.getBand().getKPointGrid()(0, 0, 0).getBandEnergy(SpinState::Up);
        std::cout << band.format() << std::endl;

        {
            std::ofstream fout("rho");
            fout << solver.getDensity();
        }
        {
            std::ofstream fout("wave");
            fout << solver.getOrbits()[0][SpinState::Up].getCoeffGrid();
        }
    }
    return 0;
}
