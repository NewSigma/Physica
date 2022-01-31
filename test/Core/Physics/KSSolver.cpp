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
#include "Physica/Utils/TestHelper.h"
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "Physica/Core/Physics/ElectronicStructure/DFT/KSSolver.h"
#include "Physica/Core/Physics/ElectronicStructure/DFT/XCProvider/LDA.h"
#include "Physica/Core/Physics/ElectronicStructure/DFT/BandGrid.h"

using namespace Physica::Core;
using ScalarType = Scalar<Double, false>;
using ComplexType = ComplexScalar<ScalarType>;

int main() {
    CrystalCell Si({5, 0, 0, 0, 5, 0, 0, 0, 5}, {0.5, 0.5, 0.5}, {14});
    ScalarType cutEnergy(0.8);
    Vector<ScalarType> data{-7.792391329, -1.041100405, -1.035201289, -1.034604466, 0.6683090416, 1.089343903, 1.092870102, 1.320171657, 1.333518296, 1.338267588, 2.048168732, 2.067794503, 2.068418852, 2.279423053, 2.296210041, 2.299139794, 2.319865956, 2.599589113, 2.607027813, 2.783839081, 3.224438445, 3.230179297, 3.239812718, 3.45520247, 3.466638718, 3.467162989, 3.613037906};
    {
        BandGrid<ScalarType, true> grid(cutEnergy, Si.reciprocal().getLattice(), 1, 1, 1, 14);
        auto solver = KSSolver<ScalarType, LDA<ScalarType, LDAType::HL, true>>(Si, cutEnergy, std::move(grid), 100, 100, 100);
        solver.solve(1E-3, 100);
        const auto& band = solver.getBand();
        Vector<ComplexType> delta = abs(band.getKPoints()[0].getEigUp().getEigenvalues() - data);
        for (size_t i = 0; i < delta.getLength(); ++i)
            if (scalarNear(delta.calc(i), ComplexType::Zero(), 1E-15))
                return 1;
    }
    {
        BandGrid<ScalarType, false> grid(cutEnergy, Si.reciprocal().getLattice(), 1, 1, 1, 14);
        auto solver = KSSolver<ScalarType, LDA<ScalarType, LDAType::HL, false>>(Si, cutEnergy, std::move(grid), 100, 100, 100);
        solver.solve(1E-3, 100);
        const auto& band = solver.getBand();
        Vector<ComplexType> delta = abs(band.getKPoints()[0].getEig().getEigenvalues() - data);
        for (size_t i = 0; i < delta.getLength(); ++i)
            if (scalarNear(delta.calc(i), ComplexType::Zero(), 1E-15))
                return 1;
    }
    return 0;
}
