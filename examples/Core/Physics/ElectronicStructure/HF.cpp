/*
 * Copyright 2021-2025 Weibo He.
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
#include "Physica/Core/Physics/ElectronicStructure/HF/RHFSolver.h"
#include "Physica/Core/Physics/ElectronicStructure/HF/GTOnG.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica;
using ScalarType = float64;

void scf_solve(double dist, ScalarType& electronEnergy, ScalarType& potentialEnergy) {
    Molecular<ScalarType> H2 = Molecular<ScalarType>(2);
    auto& atoms = H2.getAtoms();
    const Vector3D<ScalarType> pos_H1{0, 0, 0};
    const Vector3D<ScalarType> pos_H2{0, 0, dist};
    atoms[0] = pos_H1;
    atoms[1] = pos_H2;
    auto& atomicNumbers = H2.getAtomicNumbers();
    atomicNumbers[0] = 1;
    atomicNumbers[1] = 1;

    ElectronConfig config = ElectronConfig(2);
    config.setOrbitState(0, ElectronConfig::DoubleOccupacy);

    const VectorND<ScalarType> alphas{13.00773, 1.962079, 0.444529, 0.1219492, 13.00773, 1.962079, 0.444529, 0.1219492};
    RHFSolver<GaussBase<ScalarType>> solver = RHFSolver<GaussBase<ScalarType>>(H2, config, alphas.getLength());
    auto& baseSet = solver.getBaseSet();
    for (size_t i = 0; i < alphas.getLength(); ++i)
        baseSet[i] = GaussBase<ScalarType>(i < 4 ? pos_H1 : pos_H2, abs(alphas[i]), 0, 0, 0);

    if (!solver.compute(1E-5, 6)) {
        std::cerr << "[Error]: Cannot converge\n";
        exit(EXIT_FAILURE);
    }
    electronEnergy = solver.getSelfConsistentEnergy();
    potentialEnergy = H2.getNuclearRepulsionEnergy();
}
/**
 * Reference:
 * [1] J. H. Thijssen. Computational Physics[M]. London: Cambridge University Press, 2013:84
 * [2] Larsen A, Poulsen R S. Applied Hartree-Fock methods. https://projekter.aau.dk/projekter/files/213481562/report.pdf
 * 
 * TODO: Add p-orbital may improve the result, but will lead to divergence, try to enhance the DIIS procedure
 */
int main(int argc, char** argv) {
    const double from = 0.1;
    const double to = 8;
    const size_t sampleCount = 800;

    VectorND<ScalarType> r = VectorND<ScalarType>(sampleCount);
    VectorND<ScalarType> electronEnergy = VectorND<ScalarType>(sampleCount);
    VectorND<ScalarType> potentialEnergy = VectorND<ScalarType>(sampleCount);

    const double step = (to - from) / sampleCount;
    double radius = from;
    for (size_t i = 0; i < sampleCount; ++i) {
        r[i] = ScalarType(radius);
        scf_solve(radius, electronEnergy[i], potentialEnergy[i]);
        radius += step;
    }
    const VectorND<ScalarType> totalEnergy = electronEnergy + potentialEnergy;
    size_t minEnergyIndex = 0;
    for (size_t i = 1; i < sampleCount; ++i) {
        if (totalEnergy[i] > totalEnergy[minEnergyIndex])
            break;
        minEnergyIndex = i;
    }
    std::cout << "Minimum energy: " << totalEnergy[minEnergyIndex] << " At: " << r[minEnergyIndex] << std::endl;

    QApplication app(argc, argv);
    Plot* plot = new Plot(0, 8, -3, 6, 2, 2);
    auto* axisX = plot->getAxisX();
    auto* axisY = plot->getAxisY();
    axisX->setLabelFormat("%d");
    axisX->setTitleText("R/Bohr radii");
    axisY->setLabelFormat("%d");
    axisY->setTitleText("Energy/Hartree");

    plot->spline(r, electronEnergy).setName("E<sub>electron</sub>");
    plot->spline(r, potentialEnergy).setName("E<sub>nuclear</sub>");
    plot->spline(r, totalEnergy).setName("E<sub>electron</sub> + E<sub>nuclear</sub>");

    auto& chart = *plot->getChart();
    chart.setTitle("Dissociation curve of H<sub>2</sub>");
    auto* legend = chart.legend();
    legend->setAlignment(Qt::AlignTop);
    chart.setTitleFont(legend->font());
    plot->show();
    return QApplication::exec();
}
