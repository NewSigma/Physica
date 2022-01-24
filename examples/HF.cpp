/*
 * Copyright 2021 WeiBo He.
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
#include <QApplication>
#include "Physica/Core/Math/Optimization/GeneAlgorithm.h"
#include "Physica/Core/Physics/ElectronicStructure/HF/RHFSolver.h"
#include "Physica/Core/Physics/ElectronicStructure/HF/GTOnG.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Core;
using namespace Physica::Core::Physics;
using namespace Physica::Gui;

using ScalarType = Scalar<Double, false>;

void scf_solve(double dist, ScalarType& electronEnergy, ScalarType& potentialEnergy) {
    Molecular<ScalarType> H2 = Molecular<ScalarType>(2);
    auto& atoms = H2.getAtoms();
    const Vector<ScalarType, 3> pos_H1{0, 0, 0};
    const Vector<ScalarType, 3> pos_H2{0, 0, dist};
    atoms[0] = pos_H1;
    atoms[1] = pos_H2;
    auto& atomicNumbers = H2.getAtomicNumbers();
    atomicNumbers[0] = 1;
    atomicNumbers[1] = 1;

    ElectronConfig config = ElectronConfig(2);
    config.setOrbitState(0, ElectronConfig::DoubleOccupacy);
    const Vector<ScalarType, 8> alphas{13.00773, 1.962079, 0.444529, 0.1219492, 13.00773, 1.962079, 0.444529, 0.1219492};
    RHFSolver<GaussBase<ScalarType>> solver = RHFSolver<GaussBase<ScalarType>>(H2, config, alphas.getLength());
    auto& baseSet = solver.getBaseSet();
    size_t i = 0;
    baseSet[i++] = GaussBase<ScalarType>(pos_H1, abs(alphas[0]), 0, 0, 0);
    baseSet[i++] = GaussBase<ScalarType>(pos_H1, abs(alphas[1]), 0, 0, 0);
    baseSet[i++] = GaussBase<ScalarType>(pos_H1, abs(alphas[2]), 0, 0, 0);
    baseSet[i++] = GaussBase<ScalarType>(pos_H1, abs(alphas[3]), 0, 0, 0);
    baseSet[i++] = GaussBase<ScalarType>(pos_H2, abs(alphas[4]), 0, 0, 0);
    baseSet[i++] = GaussBase<ScalarType>(pos_H2, abs(alphas[5]), 0, 0, 0);
    baseSet[i++] = GaussBase<ScalarType>(pos_H2, abs(alphas[6]), 0, 0, 0);
    baseSet[i++] = GaussBase<ScalarType>(pos_H2, abs(alphas[7]), 0, 0, 0);
    if (!solver.compute(1E-5, 6)) {
        std::cerr << "[Error]: Cannot converge\n";
        exit(EXIT_FAILURE);
    }
    electronEnergy = solver.getSelfConsistentEnergy();
    potentialEnergy = H2.getNuclearRepulsionEnergy();
}
/**
 * Reference:
 * [1] Jos Thijssen. Computational Physics[M].London: Cambridge university press, 2013:84
 * [2] Larsen A, Poulsen R S. Applied Hartree-Fock methods.
 * 
 * TODO: Add p-orbital may improve the result, but will lead to divergence, try to enhance the DIIS procedure
 */
int main(int argc, char** argv) {
    const double from = 0.1;
    const double to = 8;
    const size_t sampleCount = 800;

    Vector<ScalarType> r = Vector<ScalarType>(sampleCount);
    Vector<ScalarType> electronEnergy = Vector<ScalarType>(sampleCount);
    Vector<ScalarType> potentialEnergy = Vector<ScalarType>(sampleCount);

    const double step = (to - from) / sampleCount;
    double radius = from;
    for (size_t i = 0; i < sampleCount; ++i) {
        r[i] = ScalarType(radius);
        scf_solve(radius, electronEnergy[i], potentialEnergy[i]);
        radius += step;
    }
    const Vector<ScalarType> totalEnergy = electronEnergy + potentialEnergy;
    size_t minEnergyIndex = 0;
    for (size_t i = 1; i < sampleCount; ++i) {
        if (totalEnergy[i] > totalEnergy[minEnergyIndex])
            break;
        minEnergyIndex = i;
    }
    std::cout << "Minimum energy: " << totalEnergy[minEnergyIndex] << " At: " << r[minEnergyIndex] << std::endl;

    QApplication app(argc, argv);
    Plot* plot = new Plot();

    auto& line1 = plot->spline(r, electronEnergy);
    line1.setName("Electron energy");

    auto& line2 = plot->spline(r, potentialEnergy);
    line2.setName("Nuclear repulsion energy");

    auto& line3 = plot->spline(r, totalEnergy);
    line3.setName("Effective energy");

    auto& chart = *plot->chart();
    chart.setTitle("Energy in hytrogen molecule");
    chart.legend()->setAlignment(Qt::AlignRight);
    chart.createDefaultAxes();
    chart.axes(Qt::Horizontal).first()->setTitleText("R/Bohr radii");
    chart.axes(Qt::Vertical).first()->setTitleText("Energy/Hartree");
    plot->show();
    return QApplication::exec();
}
