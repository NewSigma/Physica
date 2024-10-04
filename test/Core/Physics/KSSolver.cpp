/*
 * Copyright 2021-2023 Weibo He.
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
#include "Physica/Core/Physics/SolidState/CrystalCell.h"
#include "Physica/Core/Physics/ElectronicStructure/DFT/KSSolver.h"
#include "Physica/Core/Physics/ElectronicStructure/DFT/XCProvider/LDA.h"
#include "Physica/Core/Physics/ElectronicStructure/DFT/BandGrid.h"

using namespace Physica::Core;
using ScalarType = float64;
using ComplexType = Complex<ScalarType>;

namespace Physica {
    class Test {
    public:
        static void testCalcDensity() {
            CrystalCell<ScalarType> silicon({{5, 0, 0, 0, 5, 0, 0, 0, 5}, {0.5, 0.5, 0.5}, CrystalCell<ScalarType>::Type::Direct}, {14});
            BandGrid<ScalarType, false> bandGrid(silicon.makeRepLattice(), 1, 1, 1, 14, 8);
            KSSolver<ScalarType, LDA<ScalarType, LDAType::HL, false>> solver(silicon, 1.0, 4.5, bandGrid, 8);

            std::mt19937 gen{};
            auto& orbit = solver.getOrbits()[0][SpinState::Up];
            orbit.random_normal(gen);
            orbit.normalize();
            if (!scalarNear(orbit.calcNumElectron(), ScalarType(1), 1E-15))
                exit(EXIT_FAILURE);

            solver.calcDensity(orbit);
            auto& rSpace = solver.fft.getRSpace();
            const ScalarType numElectron = rSpace.flatten().squaredNorm();
            if (!scalarNear(numElectron, ScalarType(1), 1E-15))
                exit(EXIT_FAILURE);
        }
    };
}

void testSi() {
    constexpr double cutEnergyPsi = 0.8;
    constexpr double cutEnergyRho = 3.5;
    std::mt19937 gen{};
    CrystalCell<ScalarType> Si({{5, 0, 0, 0, 5, 0, 0, 0, 5}, {0.5, 0.5, 0.5}, CrystalCell<ScalarType>::Type::Direct}, {14});
    Vector<ScalarType> data{-7.792391329, -1.041100405, -1.035201289, -1.034604466, 0.6683090416, 1.089343903, 1.092870102, 1.320171657, 1.333518296, 1.338267588, 2.048168732, 2.067794503, 2.068418852, 2.279423053, 2.296210041, 2.299139794, 2.319865956, 2.599589113, 2.607027813, 2.783839081, 3.224438445, 3.230179297, 3.239812718, 3.45520247, 3.466638718, 3.467162989, 3.613037906};
    {
        BandGrid<ScalarType, false> grid(Si.makeRepLattice(), 1, 1, 1, 14, 8);
        auto solver = KSSolver<ScalarType, LDA<ScalarType, LDAType::HL, false>>(
                Si, cutEnergyPsi, cutEnergyRho, std::move(grid), 8);
        solver.getOrbits()[0][SpinState::Up].random_normal(gen);
        solver.getDensity().initDensity(ScalarType(Si.getElectronCount()) / Si.getVolume());
        solver.solve(1E-3, 100);
        const auto& band = solver.getBand();
        Vector<ComplexType> delta = abs(band.getKPointGrid()(0, 0, 0).getBandEnergy(SpinState::Up) - data);
        for (size_t i = 0; i < delta.getLength(); ++i)
            if (scalarNear(delta.calc(i), ComplexType(0), 1E-15))
                exit(EXIT_FAILURE);
    }
}

int main() {
    Physica::Test::testCalcDensity();
    return 0;
}
