#include <iostream>
#include <fstream>
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/Thermostat/DoubleThermo.h"
#include "Physica/Core/Physics/MD/KineticModel/PeriodicModel.h"
#include "Physica/Core/Physics/MD/ForceModel/PairModel.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Utils/Random.h"

using namespace Physica::Core;
using namespace Physica::Core::Parallel;
using ScalarType = Scalar<Double, false>;
using PosScalarType = Scalar<Double, false>;
using ThermostatType = DoubleThermo<ScalarType, PosScalarType>;
using KineticModel = PeriodicModel<ScalarType, PosScalarType, 3>;
constexpr size_t numReplica = 48;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(14);
constexpr double thermostatTime = PhyConst<AU>::secondToTime(0.1 * 1E-12);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.5;
constexpr unsigned int numMolecular = 108;
constexpr double pair_cutoff = 15;
constexpr double molarVolume = 25.6;
constexpr double mass = PhyConst<AU>::atomMass(1) * 2;

ScalarType pot_functor(ScalarType r) {
    constexpr double alpha = 1.713;
    constexpr double beta = 1.5671;
    constexpr double gamma = 0.00993;
    constexpr double cutoff = 8.32;
    constexpr double c6 = 12.14;
    constexpr double c8 = 215.2;
    constexpr double c9 = 143.1;
    constexpr double c10 = 4813.9;

    const ScalarType r2 = square(r);
    ScalarType result = exp(-r2 * gamma - r * beta + alpha);
    const ScalarType rep_r = reciprocal(r);
    const ScalarType rep_r2 = square(rep_r);
    const ScalarType rep_r4 = square(rep_r2);
    const ScalarType rep_r6 = rep_r4 * rep_r2;
    const ScalarType rep_r8 = square(rep_r4);
    const ScalarType rep_r9 = rep_r8 * rep_r;
    const ScalarType rep_r10 = rep_r6 * rep_r4;
    const ScalarType g = rep_r6 * c6 + rep_r8 * c8 - rep_r9 * c9 + rep_r10 * c10;

    if (r < cutoff) {
        const ScalarType f_cutoff = exp(-square(rep_r * cutoff - 1));
        result -= g * f_cutoff;
    }
    else
        result -= g;
    return result;
}

ScalarType force(ScalarType r) {
    constexpr double alpha = 1.713;
    constexpr double beta = 1.5671;
    constexpr double gamma = 0.00993;
    constexpr double cutoff = 8.32;
    constexpr double c6 = 12.14;
    constexpr double c8 = 215.2;
    constexpr double c9 = 143.1;
    constexpr double c10 = 4813.9;

    const ScalarType r2 = square(r);
    ScalarType result = exp(-r2 * gamma - r * beta + alpha) * (r * (gamma * 2) + beta);
    const ScalarType rep_r = reciprocal(r);
    const ScalarType rep_r2 = square(rep_r);
    const ScalarType rep_r4 = square(rep_r2);
    const ScalarType rep_r6 = rep_r4 * rep_r2;
    const ScalarType rep_r7 = rep_r6 * rep_r;
    const ScalarType rep_r9 = rep_r7 * rep_r2;
    const ScalarType rep_r10 = rep_r6 * rep_r4;
    const ScalarType rep_r11 = rep_r7 * rep_r4;

    const ScalarType d_g = rep_r10 * (9 * c9) - rep_r11 * (10 * c10) - rep_r9 * (8 * c8) - rep_r7 * (6 * c6);
    if (r < cutoff) {
        const ScalarType rep_r8 = square(rep_r4);
        const ScalarType g = rep_r6 * c6 + rep_r8 * c8 - rep_r9 * c9 + rep_r10 * c10;
        const ScalarType f_cutoff = exp(-square(rep_r * cutoff - 1));
        result += (d_g + g * rep_r * (rep_r2 * (2 * cutoff * cutoff) - rep_r * (2 * cutoff))) * f_cutoff;
    }
    else {
        result += d_g;
    }
    return result;
}

template<class RandomGenerator>
RPMD<ScalarType, PosScalarType> makeSystem(RandomGenerator& gen) {
    using MDCellType = typename RPMD<ScalarType, PosScalarType>::MDCellType;
    typename MDCellType::LatticeMatrix lattice = MDCellType::LatticeMatrix::unitMatrix(3);
    typename MDCellType::PositionMatrix pos(numMolecular, 3);
    std::uniform_real_distribution dist{};
    for (auto& elem : pos)
        elem = dist(gen);
    typename MDCellType::MassVector massVec(numMolecular, mass);
    MDCellType cell(std::move(lattice), std::move(pos), std::move(massVec));

    const double factor = (std::cbrt(numMolecular * molarVolume / PhyConst<SI>::avogadroNa) / 100) / PhyConst<SI>::bohrRadius;
    cell.scale(factor);

    return RPMD<ScalarType, PosScalarType>(std::move(cell), numReplica, numReplica, temperatureT, timeStep);
}
/**
 * Reference:
 * [1] Miller TF, Manolopoulos DE. 2005. Quantum diffusion in liquid para-hydrogen from ring polymer molecular dynamics. J. Chem. Phys. 122:184503
 */
int main() {
    constexpr size_t CorrStep = 3000;
    constexpr double factor = 1.0 / (numMolecular * mass * mass) * (PhyConst<AU>::bohrToAngstorm(1) * PhyConst<AU>::bohrToAngstorm(1)) / (PhyConst<AU>::timeToSecond(1) * 1E12 * PhyConst<AU>::timeToSecond(1) * 1E12);

    ThreadPool::initThreadPool(4);
    ThreadPool& pool = ThreadPool::getInstance();
    {
        std::mt19937 gen(1082247429173841685);
        RPMD<ScalarType, PosScalarType> rpmd = makeSystem(gen);
        auto& ringPolymer = rpmd.getRingPolymer();
        const ThermostatType thermo(temperatureT, thermostatTime);
        rpmd.initMomentum(gen);

        PairModel<ScalarType, PosScalarType, decltype(&force)> pair(ScalarType(pair_cutoff), force, pot_functor);
        rpmd.updateForce<decltype(pair), ThreadExecutor>(pair);
        KineticModel kineticModel(temperatureT, numReplica);

        Vector<ScalarType> corr(CorrStep, 0);
        Vector<ScalarType> corr_var(CorrStep, 0);
        Vector<ScalarType> temp(CorrStep);

        for (unsigned int i = 0; i < 6; ++i) {
            temp = ScalarType(0);
            for (unsigned int j = 0; j < 100; ++j) {
                rpmd.nvt_step_for<ThermostatType, decltype(gen), KineticModel, decltype(pair), ThreadExecutor>(PhyConst<AU>::secondToTime(2 * 1E-12), thermo, gen, kineticModel, pair);
                const auto p0 = ringPolymer.makeCentroidMomentum();
                for (unsigned int k = 0; k < CorrStep; ++k) {
                    toNextMean(temp[k], j, hadamard(ringPolymer.makeCentroidMomentum(), p0).sum() * factor);
                    rpmd.nvt_step<ThermostatType, decltype(gen), KineticModel, decltype(pair), ThreadExecutor>(thermo, gen, kineticModel, pair);
                }
            }

            for (unsigned int k = 0; k < CorrStep; ++k)
                toNextVariance(corr_var[k], corr[k], i, temp[k]);
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
    pool.shouldExit();
    ThreadPool::deInitThreadPool();
    return 0;
}
