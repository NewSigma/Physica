#include <iostream>
#include <fstream>
#include "Physica/Core/Math/Random/RandomSeed.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/Thermostat/DoubleThermo.h"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Core/Physics/MD/ForceModel/SilveraGoldman.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"

using namespace Physica::Core;
using ScalarType = Scalar<Double>;
using ThermostatType = DoubleThermo<ScalarType>;
using KineticModel = FreeModel<ScalarType, 3, Dynamic, RPMDIntegrator::Exact>;
using ForceModel = SilveraGoldman<ScalarType>;
using RandomPoolType = RandomPool<std::mt19937>;
constexpr size_t numReplica = 48;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(14);
constexpr double thermostatTime = PhyConst<AU>::secondToTime(0.1 * 1E-12);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.5;
constexpr unsigned int numMolecular = 108;
constexpr double pair_cutoff = 15;
constexpr double molarVolume = 25.6;
constexpr double mass = PhyConst<AU>::atomMass(1) * 2;
constexpr size_t CorrStep = 3000;

template<class RandomGenerator>
RPMD<ScalarType> makeSystem(RandomGenerator& gen) {
    using MDCellType = typename RPMD<ScalarType>::MDCellType;
    typename MDCellType::LatticeMatrix lattice = MDCellType::LatticeMatrix::unitMatrix(3);
    typename MDCellType::PositionMatrix pos(numMolecular, 3);
    std::uniform_real_distribution dist{};
    for (auto& elem : pos)
        elem = dist(gen);
    typename MDCellType::MassVector massVec(numMolecular, mass);
    MDCellType cell(std::move(lattice), std::move(pos), std::move(massVec));

    const double factor = (std::cbrt(numMolecular * molarVolume / PhyConst<SI>::avogadroNa) / 100) / PhyConst<SI>::bohrRadius;
    cell.scale(factor);

    return RPMD<ScalarType>(std::move(cell), numReplica, numReplica, temperatureT, timeStep);
}
/**
 * Reference:
 * [1] Miller TF, Manolopoulos DE. 2005. Quantum diffusion in liquid para-hydrogen from ring polymer molecular dynamics. J. Chem. Phys. 122:184503
 */
int main() {
    ThreadPool::numThreadRequired = 4;
    Vector<ScalarType> corr(CorrStep, 0);
    Vector<ScalarType> corr_var(CorrStep, 0);
    {
        constexpr double factor = 1.0 / (numMolecular * mass * mass) * (PhyConst<AU>::bohrToAngstorm(1) * PhyConst<AU>::bohrToAngstorm(1)) / (PhyConst<AU>::timeToSecond(1) * 1E12 * PhyConst<AU>::timeToSecond(1) * 1E12);
        auto& gen = RandomPoolType::getGen();
        RPMD<ScalarType> rpmd = makeSystem(gen);
        rpmd.initMomentum(gen);

        ForceModel forceModel(pair_cutoff);
        KineticModel kineticModel(temperatureT, numReplica);
        const ThermostatType thermo(temperatureT, thermostatTime);

        auto& ringPolymer = rpmd.getRingPolymer();
        Vector<ScalarType> temp(CorrStep);
        for (unsigned int i = 0; i < 6; ++i) {
            temp = ScalarType(0);
            for (unsigned int j = 0; j < 100; ++j) {
                rpmd.nvt_step_for<ThermostatType, RandomPoolType, KineticModel, ForceModel, ThreadExecutor>(PhyConst<AU>::secondToTime(2 * 1E-12), thermo, kineticModel, forceModel);
                const auto p0 = ringPolymer.makeCentroidMomentum();
                for (unsigned int k = 0; k < CorrStep; ++k) {
                    toNextMean(temp[k], j, hadamard(ringPolymer.makeCentroidMomentum(), p0).sum() * factor);
                    rpmd.nvt_step<ThermostatType, RandomPoolType, KineticModel, ForceModel, ThreadExecutor>(thermo, kineticModel, forceModel);
                }
            }
            toNextVariance(corr_var, corr, i, temp);
        }

        {
            std::ofstream fout("corr");
            fout << corr;
        }
        {
            std::ofstream fout("corr_var");
            fout << corr_var;
        }
    }
    ThreadPool::getInstance().shouldExit();
    return 0;
}
