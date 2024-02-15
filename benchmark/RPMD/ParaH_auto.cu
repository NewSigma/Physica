/*
 * Copyright 2023-2024 WeiBo He.
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
#include <gperftools/profiler.h>
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/Thermostat/TRPMDThermo.h"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Core/Physics/MD/ForceModel/SilveraGoldman.cuh"
#include "Physica/Core/Parallel/Executor/AutoExecutor.h"
#include "Physica/Core/Parallel/Executor/CudaExecutor.cuh"
#include "Physica/Core/Math/Random/RandomPool.h"
#include "Physica/Utils/BenchmarkHelper.h"

using namespace Physica::Core;
using ScalarType = Scalar<Float>;
using MDType = RPMD<ScalarType, 3, Physica::Utils::Dynamic, Physica::Utils::PageLockedAllocator<ScalarType>>;
using MDCellType = typename MDType::MDCellType;
using KineticModel = FreeModel<ScalarType, 3, Dynamic, RPMDIntegrator::Exact>;
using RandomPoolType = RandomPool<std::mt19937, 10000>;
constexpr size_t numReplica = 24;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(25);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.5;
constexpr double pair_cutoff = 15;
constexpr double molarVolume = 31.7;
constexpr double mass = PhyConst<AU>::atomMass(1) * 2;

template<bool IsSmallCell>
class MixedModel {
    using DeviceModelType = device_obj<SilveraGoldman<ScalarType, true, IsSmallCell>>;
    SilveraGoldman<ScalarType, true, false> hostModel;
    Physica::Utils::Array<DeviceModelType> deviceModels;
public:
    MixedModel(size_t numParticle, ScalarType cutoff, size_t numCudaThread)
        : hostModel(cutoff)
        , deviceModels(numCudaThread, numParticle, cutoff) {
        assert(numCudaThread < ThreadPool::getInstance().getNumThreads());
    }

    template<class Executor>
    [[nodiscard]] Vector<ScalarType> force(const MDCellType& cell) {
        Vector<ScalarType> result(cell.getDOF());
        forceAsync<Vector<ScalarType>, Executor>(cell, result);
        StreamPool::getStream().wait();
        return result;
    }

    template<class VectorType, class Executor>
    void forceAsync(const MDCellType& cell, ContinuousVector<VectorType>& result) {
        static_assert(Internal::Traits<Executor>::isCudaEnabled, "[Error]: Invalid executor");
        const auto threadId = ThreadPool::getThreadInfo().id;
        const bool useCPU = ThreadPool::isMainThread() || threadId >= getNumCudaThread();
        if (useCPU)
            result = hostModel.template force<ThreadExecutor>(cell);
        else {
            deviceModels[threadId].template forceAsync<VectorType, CudaExecutor>(cell, result);
        }
    }

    template<class Executor>
    [[nodiscard]] Vector<ScalarType> force_short(const MDCellType& cell) { return force<Executor>(cell); }
    template<class Executor>
    [[nodiscard]] Vector<ScalarType> force_long(const MDCellType& cell) const { return Vector<ScalarType>(cell.getNumParticle() * 3, 0); }
    /* Getters */
    [[nodiscard]] size_t getNumCudaThread() const noexcept { return deviceModels.getLength(); }
private:
    /* Getters */
    [[nodiscard]] bool shouldUseCPU() const noexcept {
        return ThreadPool::isMainThread() || ThreadPool::getThreadInfo().id >= getNumCudaThread();
    }
};

namespace Physica::Core::Internal {
    template<bool IsSmallCell>
    class Traits<MixedModel<IsSmallCell>> {
    public:
        constexpr static bool IsPeriodBoundary = true;
        constexpr static bool IsContractable = false;
    };
}

template<class RandomGenerator>
MDType makeSystem(size_t numMolecular, RandomGenerator& gen) {
    typename MDCellType::LatticeMatrix lattice = MDCellType::LatticeMatrix::unitMatrix(3);
    typename MDCellType::PositionMatrix pos(numMolecular, 3);
    std::uniform_real_distribution dist{};
    for (auto& elem : pos)
        elem = dist(gen);
    typename MDCellType::MassVector massVec(numMolecular, mass);
    MDCellType cell(std::move(lattice), std::move(pos), std::move(massVec));

    const double factor = (std::cbrt(numMolecular * molarVolume / PhyConst<SI>::avogadroNa) / 100) / PhyConst<SI>::bohrRadius;
    cell.scale(factor);

    return MDType(std::move(cell), numReplica, numReplica, temperatureT, timeStep);
}
/**
 * Reference:
 * [1] Miller TF, Manolopoulos DE. 2005. Quantum diffusion in liquid para-hydrogen from ring polymer molecular dynamics. J. Chem. Phys. 122:184503
 */
int main() {
    ThreadPool::numThreadRequired = 8;
    auto& gen = RandomPoolType::getInstance().getGen();
    KineticModel kineticModel(temperatureT, numReplica);
    {
        using ForceModel = MixedModel<true>;
        constexpr size_t numMolecular = 108;
        MDType rpmd = makeSystem(numMolecular, gen);
        rpmd.initMomentum(gen);
        ForceModel forceModel(numMolecular, pair_cutoff, 4);
        auto timeuse = Physica::Utils::Benchmark::run([&]() {
            rpmd.nve_step_for<KineticModel, ForceModel, AutoExecutor>(PhyConst<AU>::secondToTime(2 * 1E-13), kineticModel, forceModel);
        }, 8, 20);
        std::cout << "108 atom time use: " << timeuse.first << '(' << timeuse.second << ")\n";
    }
    {
        using ForceModel = MixedModel<true>;
        constexpr size_t numMolecular = 256;
        MDType rpmd = makeSystem(numMolecular, gen);
        rpmd.initMomentum(gen);
        ForceModel forceModel(numMolecular, pair_cutoff, 4);
        auto timeuse = Physica::Utils::Benchmark::run([&]() {
            rpmd.nve_step_for<KineticModel, ForceModel, AutoExecutor>(PhyConst<AU>::secondToTime(2 * 1E-13), kineticModel, forceModel);
        }, 8, 20);
        std::cout << "256 atom time use: " << timeuse.first << '(' << timeuse.second << ")\n";
    }
    {
        using ForceModel = MixedModel<true>;
        constexpr size_t numMolecular = 500;
        MDType rpmd = makeSystem(numMolecular, gen);
        rpmd.initMomentum(gen);
        ForceModel forceModel(numMolecular, pair_cutoff, 4);
        auto timeuse = Physica::Utils::Benchmark::run([&]() {
            rpmd.nve_step_for<KineticModel, ForceModel, AutoExecutor>(PhyConst<AU>::secondToTime(1 * 1E-13), kineticModel, forceModel);
        }, 8, 20);
        std::cout << "500 atom time use: " << timeuse.first << '(' << timeuse.second << ")\n";
    }
    {
        using ForceModel = MixedModel<false>;
        constexpr size_t numMolecular = 864;
        MDType rpmd = makeSystem(numMolecular, gen);
        rpmd.initMomentum(gen);
        ForceModel forceModel(numMolecular, pair_cutoff, 4);
        auto timeuse = Physica::Utils::Benchmark::run([&]() {
            rpmd.nve_step_for<KineticModel, ForceModel, AutoExecutor>(PhyConst<AU>::secondToTime(5 * 1E-14), kineticModel, forceModel);
        }, 8, 20);
        std::cout << "864 atom time use: " << timeuse.first << '(' << timeuse.second << ")\n";
    }
    return 0;
}
