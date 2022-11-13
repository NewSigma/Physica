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
#include <fstream>
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "Physica/Core/Physics/ElectronicStructure/DFT/KSSolver.h"
#include "Physica/Core/Physics/ElectronicStructure/DFT/XCProvider/LDA.h"
#include "Physica/Core/Physics/ElectronicStructure/DFT/BandGrid.h"
#include "Physica/Utils/Cycler.h"
#include "Physica/Utils/Random.h"

using namespace Physica::Core;
using ScalarType = Scalar<Double, false>;
using ComplexType = ComplexScalar<ScalarType>;

int main() {
    using namespace Physica::Utils;
    Cycler::init();
    CrystalCell Si({5, 0, 0, 0, 5, 0, 0, 0, 5}, {0.5, 0.5, 0.5}, {2}, CrystalCell::Type::Direct);
    ScalarType cutEnergy(0.8);

    std::mt19937::result_type seed;
    Physica::Utils::Random::rdrand(seed);
    std::mt19937 gen(seed);
    {
        BandGrid<ScalarType, false> grid(Si.reciprocal().getLattice(), 1, 1, 1, 2, 2);
        using SolverType = KSSolver<ScalarType, LDA<ScalarType, LDAType::HL, false>>;
        SolverType solver = SolverType(Si, cutEnergy, 2, std::move(grid), gen);
/*
        using DensityType = typename SolverType::DensityType;
        DensityType rho(solver.getFFTDim());
        DensityType last_rho{};
        std::ifstream fin("rho");
        fin >> last_rho;
        rho.fit(last_rho);
        solver.setDensity(std::move(rho));
*/
        const auto from = Cycler::tic();
        solver.solve(1E-4, 100);
        const auto to = Cycler::toc();
        std::cout << "Time use: " << Cycler::toSeconds(to - from) << '\n';
        const auto& band = solver.getBand().getKPointGrid()(0, 0, 0).getBandEnergy(SpinState::Up);
        std::cout << band.format() << std::endl;

        //std::ofstream fout("rho");
        //fout << solver.getDensityGrid();
    }
    return 0;
}
