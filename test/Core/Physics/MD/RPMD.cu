/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/Thermostat/DoubleThermo.h"
#include "Physica/Core/Physics/MD/ForceModel/SilveraGoldman.cuh"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Core/Parallel/Executor/CUDAExecutor.cuh"
#include "Test.h"

using namespace Physica;
using T = float32;
using ForceModel = device_obj<SilveraGoldman<T, true>>;
using KineticModel = FreeModel<T, 3, Dynamic, RPMDIntegrator::Exact>;
using ThermoType = DoubleThermo<KineticModel>;
using MDType = RPMD<T, 3, Physica::Dynamic, PageLockedAllocator<T>>;
using RandomSource = Random<MT19937, 3438603950906262892>;
constexpr size_t numReplica = 24;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(25);
constexpr double thermostatTime = PhyConst<AU>::secondToTime(100 * 1E-15);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.5;
constexpr unsigned int numMolecular = 108;
constexpr double pair_cutoff = 15;
constexpr double molarVolume = 31.7;
constexpr double mass = PhyConst<AU>::atomMass(1) * 2;

namespace {
    MDType makeSystem() {
        using MDCellType = MDType::MDCellType;
        MDCellType::LatticeMatrix lattice = MDCellType::LatticeMatrix::identity(3);
        auto pos = MDCellType::PositionMatrix::random_uniform<RandomSource>(numMolecular, 3);
        MDCellType::MassVector massVec(numMolecular, mass);
        MDCellType cell(std::move(lattice), std::move(pos), std::move(massVec));

        const double factor = (std::cbrt(numMolecular * molarVolume / PhyConst<SI>::avogadroNa) / 100) / PhyConst<SI>::bohrRadius;
        cell.scale(factor);

        return MDType(std::move(cell), numReplica, numReplica, temperatureT, timeStep);
    }
    /**
    * Reference:
    * [1] J. Chem. Phys. 122, 184503 (2005); https://doi.org/10.1063/1.1893956
    */
    void testMDRun() {
        T mean = 0;
        T var = 0;
        {
            const ThermoType thermo(temperatureT, thermostatTime);
            KineticModel kineticModel(temperatureT, numReplica);
            ForceModel forceModel(numMolecular, pair_cutoff);
            auto rpmd = makeSystem();
            rpmd.initMomentum<KineticModel, RandomSource>();

            for (unsigned int i = 0; i < 6; ++i) {
                T temp = 0;
                rpmd.nvt_step_for<RandomSource, GPU>(
                    PhyConst<AU>::secondToTime(2 * 1E-12),
                    thermo,
                    kineticModel,
                    forceModel);

                for (unsigned int j = 0; j < 100; ++j) {
                    rpmd.nvt_step<RandomSource, GPU>(thermo, kineticModel, forceModel);
                   temp.toNextMean(j, rpmd.calcKinetic<KineticModel>());
                }
                var.toNextVariance(mean, i, temp);
            }
        }
        constexpr double answer = 61.8;
        const T energyPerMol = PhyConst<AU>::temperatureToK(double(mean) / numMolecular);
        const T deviation = PhyConst<AU>::temperatureToK(std::sqrt(double(var))) / numMolecular;
        expect(abs(energyPerMol - answer) < deviation * 2.0);
    }
}

int main() {
    testMDRun();
    return 0;
}
