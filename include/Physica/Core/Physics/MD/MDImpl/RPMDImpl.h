/*
 * Copyright 2022-2023 WeiBo He.
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
#pragma once

#include "Physica/Core/Physics/MD/ForceModel/EmptyForceModel.h"
#include "Physica/Core/Math/Random/RandomPool.h"

namespace Physica::Core {
    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::RPMD(
            MDCellType cell_,
            size_t numReplica,
            size_t numContract,
            ScalarType temperatureT_,
            ScalarType timeStep_)
            : cell(std::move(cell_))
            , fftContract(numContract, 1, PlanFlag::Estimate)
            , timeStep(std::move(timeStep_)) {
        assert(0 < numContract && numContract <= numReplica);
        assert(NumReplica == Dynamic || NumReplica == numReplica);
        ringPolymer = RingPolymerType(cell, numReplica);

        const size_t dof = getDOF();
        forceBuffer.resize(dof, numReplica);
        if (isContractEnabled()) {
            posContract.resize(dof, numContract);
            forceContract.resize(dof, numContract);
        }

        setTemperature(temperatureT_);
        checkParam();
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>&
    RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::operator=(
            RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator> obj) noexcept {
        swap(obj);
        return *this;
    }
    /**
     * Contract method to improve performance in [1]
     * Reference:
     * [1] T. E. Markland, D. E. Manolopoulos. An efficient ring polymer contraction scheme for imaginary time path integral simulations[J]. J. Chem. Phys. 129, 024105 (2008)
     */
    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class ForceModel, class Executor>
    void RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::updateForce(ForceModel& model) {
        constexpr bool isCudaEnabled = Internal::Traits<Executor>::isCudaEnabled;
        static_assert(!isCudaEnabled || std::allocator_traits<ForceMatrixAllocator>::isPageLocked
                , "[Error]: Allocator is not page locked, performance will decrease");
        if (isContractEnabled()) {
            auto kernel_short = [&](unsigned int thread) {
                const auto range = Executor::splitJob(getNumReplica(), Executor::getNumThread(), thread);
                for (size_t replica = range.first; replica < range.second; ++replica) {
                    MDCellType cell = phaseToCell(replica);
                    cell.normalize();
                    auto saveTo = forceBuffer.col(replica);
                    saveTo = model.template force_short<Executor>(std::move(cell));
                }
            };
            auto future_short = Executor::parallel_for(kernel_short, Executor::getNumThread());

            contract();
            auto kernel_long = [&](unsigned int thread) {
                const auto range = Executor::splitJob(getNumReplica(), Executor::getNumThread(), thread);
                for (size_t contract = range.first; contract < range.second; ++contract) {
                    MDCellType cell = contractToCell(contract);
                    cell.normalize();
                    auto saveTo = forceContract.col(contract);
                    saveTo = model.template force_long<Executor>(std::move(cell));
                }
            };
            auto future_long = Executor::parallel_for(kernel_long, Executor::getNumThread());

            Executor::auto_wait(future_long);
            decontract();
            Executor::auto_wait(future_short);
        }
        else {
            auto kernel = [this, &model](unsigned int replica) {
                MDCellType cell = phaseToCell(replica);
                cell.normalize();
                auto saveTo = forceBuffer.col(replica);
                using VectorType = typename decltype(saveTo)::VectorBase;
                model.template forceAsync<VectorType, Executor, false>(std::move(cell), saveTo);
            };
            auto future = Executor::parallel_for(kernel, getNumReplica());
            Executor::auto_wait(future);
            Executor::wait();
        }
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel,
             class ForceModel,
             class Executor>
    void RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::nve_step(KineticModel& kineticModel, ForceModel& forceModel) {
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        if (isFreeModel) {
            kineticModel.nve_step(ringPolymer, timeStep);
        }
        else {
            forceStep(timeStep * PlainScalar(0.5));
            kineticModel.nve_step(ringPolymer, timeStep);
            updateForce<ForceModel, Executor>(forceModel);
            forceStep(timeStep * PlainScalar(0.5));
        }
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel,
             class ForceModel,
             class Executor>
    void RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::nve_step_for(
            ScalarType duration, KineticModel& kineticModel, ForceModel& forceModel) {
        const uint64_t step = Base::durationToStep(duration, timeStep);
        for (uint64_t _ = 0; _ < step; ++_)
            nve_step<KineticModel, ForceModel, Executor>(kineticModel, forceModel);
    }
    /**
     * BAOAB integrator as introduced in [1]
     * Reference:
     * [1] Liu J, Li D, Liu X. A simple and accurate algorithm for path integral molecular dynamics with the Langevin thermostat[J]. J. Chem. Phys, 2016, 145(2):1291-1301.
     */
    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class Thermostat,
             class RandomPoolType,
             class KineticModel,
             class ForceModel,
             class Executor>
    void RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::nvt_step(
            const Thermostat& thermostat,
            KineticModel& kineticModel,
            ForceModel& forceModel) {
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        constexpr bool isSeedFixed = RandomPoolType::isSeedFixed();
        using ThermoStepExecutor = typename std::conditional<isSeedFixed, SequentialExecutor, Executor>::type;
        if (isFreeModel) {
            kineticModel.nve_step(ringPolymer, timeStep * 0.5);
            thermostat.template step<RandomPoolType, ThermoStepExecutor>(ringPolymer, timeStep);
            kineticModel.nve_step(ringPolymer, timeStep * 0.5);
        }
        else {
            forceStep(timeStep * 0.5);
            kineticModel.nve_step(ringPolymer, timeStep * 0.5);
            thermostat.template step<RandomPoolType, ThermoStepExecutor>(ringPolymer, timeStep);
            kineticModel.nve_step(ringPolymer, timeStep * 0.5);
            updateForce<ForceModel, Executor>(forceModel);
            forceStep(timeStep * 0.5);
        }
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class Thermostat,
             class RandomPoolType,
             class KineticModel,
             class ForceModel,
             class Executor>
    void RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::nvt_step_for(
            ScalarType duration,
            const Thermostat& thermostat,
            KineticModel& kineticModel,
            ForceModel& forceModel) {
        const uint64_t step = Base::durationToStep(duration, timeStep);
        for (uint64_t _ = 0; _ < step; ++_)
            nvt_step<Thermostat, RandomPoolType, KineticModel, ForceModel, Executor>(thermostat, kineticModel, forceModel);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class Thermostat,
             class RandomGenerator,
             class Barostat,
             class KineticModel,
             class ForceModel,
             class Executor>
    void RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::npt_step(
            const Thermostat& thermostat,
            RandomGenerator& gen,
            Barostat& barostat,
            KineticModel& kineticModel,
            ForceModel& forceModel) {
        barostat.forceStep(*this, timeStep * 0.5);
        kineticModel.npt_step(ringPolymer, cell, barostat, timeStep * 0.5);
        thermostat.step(ringPolymer, gen, timeStep);
        kineticModel.npt_step(ringPolymer, cell, barostat, timeStep * 0.5);
        updateForce<ForceModel, Executor>(forceModel);
        barostat.forceStep(*this, timeStep * 0.5);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class Thermostat, class RandomGenerator, class Barostat, class KineticModel, class ForceModel, class Executor>
    void RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::npt_step_for(
            ScalarType duration,
            const Thermostat& thermostat,
            RandomGenerator& gen,
            Barostat& barostat,
            KineticModel& kineticModel,
            ForceModel& forceModel) {
        const uint64_t step = Base::durationToStep(duration, timeStep);
        for (uint64_t _ = 0; _ < step; ++_)
            npt_step<Thermostat, RandomGenerator, Barostat, KineticModel, ForceModel, Executor>(thermostat, gen, barostat, kineticModel, forceModel);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class RandomGenerator>
    void RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::initMomentum(RandomGenerator& gen) {
        return ringPolymer.initMomentum(temperatureT, gen);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::scaleVelocity() {
        ringPolymer.scaleVelocity(temperatureT);
    }
    /**
     * Carrying out this function every several steps may stable the simulation.
     */
    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::normalizeCentroid() {
        PositionMatrix centroid = ringPolymer.makeCentroidPos();
        cell.toDirect(centroid);
        size_t index = getDOF();
        for (const auto elem : centroid) {
            const size_t component = index % Dim;
            const size_t atom_start = index - component;
            const int integer = float(elem);
            const Vector<PosScalarType, Dim> delta = PosScalarType(integer - elem.isNegative()) * cell.getLattice().row(component).asVector();
            for (size_t i = 0; i < Dim; ++i) {
                auto row = getPhaseMatrix().row(atom_start + i);
                row -= delta[i];
            }
            ++index;
        }
        assert(checkCentroid());
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    typename RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::MDCellType
    RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::phaseToCell(size_t replica) const {
        return MDCellType(cell.getLattice(), ringPolymer.makeBeadPos(replica), cell.getMassVec());
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    typename RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::MDCellType
    RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::contractToCell(size_t contract) const {
        PositionMatrix pos(getNumParticle(), Dim);
        auto phase = posContract.col(contract);
        size_t index = 0;
        for (auto& elem : pos) {
            elem = PosScalarType(phase[index]);
            ++index;
        }
        return MDCellType(cell.getLattice(), std::move(pos), cell.getMassVec());
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    typename RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::MDCellType
    RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::makeAverageCell() const {
        return MDCellType(getLattice(), ringPolymer.makeCentroidPos(), cell.getMassVec());
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::swap(RPMD& obj) noexcept {
        cell.swap(obj.cell);
        ringPolymer.swap(obj.ringPolymer);
        forceBuffer.swap(obj.forceBuffer);

        fftContract.swap(obj.fftContract);
        posContract.swap(obj.posContract);
        forceContract.swap(obj.forceContract);

        temperatureT.swap(obj.temperatureT);
        timeStep.swap(obj.timeStep);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class ForceModel, class Executor>
    ScalarType RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::calcPotential(const ForceModel& model) const {
        Vector<ScalarType> temp(getNumReplica());
        auto kernel = [this, model, &temp](unsigned int replica) {
            temp[replica] = model.potentialEnergy(phaseToCell(replica));
        };
        Executor::parallel_for(kernel, getNumReplica(), Executor::getNumThread()).wait();
        return mean(temp);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class ForceModel>
    typename RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::LatticeMatrix
    RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::makeStress(const ForceModel& model) const {
        const size_t dof = getDOF();
        LatticeMatrix stress(Dim, Dim, 0);
        for (size_t replica = 0; replica < 1; ++replica) {
            const auto col = getPhaseMatrix().col(replica);
            const auto momentum = col.head(dof);
            const auto pos = col.tail(dof);
            const auto force = forceBuffer.col(replica);
            for (size_t i = 0; i < getNumParticle(); ++i) {
                const size_t from = i * getDim();
                const size_t to = from + getDim();
                const ScalarType factor = reciprocal(getMassVec()[i]);
                const auto momentum1 = momentum.segment(from, to);
                const auto pos1 = pos.segment(from, to);
                const auto force1 = force.segment(from, to);
                stress += (factor * momentum1) * momentum1.transpose();
            }
        }
        stress *= reciprocal(getVolume());
        return stress + model.virial(makeAverageCell());
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class ForceModel>
    ScalarType RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::calcClassicalPotentialEnergy(const ForceModel& model) const {
        ScalarType result = 0;
        for (size_t i = 0; i < getNumReplica(); ++i)
            result += model.potentialEnergy(phaseToCell(i));
        return result;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class ForceModel>
    ScalarType RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::calcClassicalInternalEnergy(const ForceModel& model) const {
        return ringPolymer.calcClassicalKinetic() + calcClassicalPotentialEnergy<ForceModel>(model) + calcClassicalElastic();
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    ScalarType RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::calcClassicalElastic() const {
        const size_t dof = getDOF();
        auto pos = getPhaseMatrix().bottomRows(dof);
        const ScalarType omegaW = ringPolymer.calcOmegaW(temperatureT);
        ScalarType result = 0;
        for (size_t i = 0; i < dof; ++i) {
            const ScalarType mass = cell.getMass(i / Dim);
            for (size_t j = 0; j < getNumReplica(); ++j)
                result += mass * square(omegaW * (pos(i, j) - pos(i, (j + 1) % getNumReplica()))) * 0.5;
        }
        return result;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    ScalarType RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::calcKinetic() const {
        const ScalarType repBeta = ringPolymer.calcRepBeta(ringPolymer.calcTemperature());
        const size_t dof = getDOF();
        const auto centroidPos = ringPolymer.makeCentroidPos();

        ScalarType kinetic = repBeta * dof;
        for (size_t replica = 0; replica < getNumReplica(); ++replica) {
            auto phase = getPhaseMatrix().col(replica);
            auto pos = phase.tail(dof);
            kinetic += (centroidPos.flatten() - pos) * forceBuffer.col(replica);
        }
        kinetic /= ScalarType(getNumReplica() * 2);
        return kinetic;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::setTemperature(ScalarType temperature) {
        assert(!temperature.isNegative());
        temperatureT = temperature;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::toContractBeadRepr(size_t posID) {
        assert(posID < getDOF());
        auto row = ringPolymer.getBuffer().row(1);
        auto head = row.head(fftContract.getKSpaceSize());
        fftContract.invTransform(head);
        auto pos = posContract.row(posID);
        pos = fftContract.getRSpace() * (ScalarType(getNumContract()) / ScalarType(getNumReplica()));
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::forceToNormRepr(size_t posID) {
        assert(posID < getDOF());
        fftContract.transform(forceContract.row(posID));
        auto row = ringPolymer.getBuffer().row(0);
        row.asVector() = ScalarType(0);
        auto head = row.head(fftContract.getKSpaceSize());
        head = fftContract.getKSpace();
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::forceToBeadRepr(size_t posID) {
        assert(posID < getDOF());
        auto& buffer = ringPolymer.getBuffer();
        auto& fft = ringPolymer.getCanonicalFFT();
        fft.invTransform(buffer.row(0));
        auto f = forceBuffer.row(posID);
        f.asVector() += fft.getRSpace() * (ScalarType(getNumReplica()) / ScalarType(getNumContract()));
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::contract() {
        const size_t dof = getDOF();
        for (size_t i = 0; i < dof; ++i) {
            ringPolymer.toNormalRepr(i);
            toContractBeadRepr(i);
        }
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::decontract() {
        const size_t dof = getDOF();
        for (size_t i = 0; i < dof; ++i) {
            forceToNormRepr(i);
            forceToBeadRepr(i);
        }
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::forceStep(ScalarType deltaT) {
        const size_t dof = getDOF();
        for (size_t replica = 0; replica < getNumReplica(); ++replica) {
            auto col = getPhaseMatrix().col(replica);
            auto momentum = col.head(dof);
            momentum += forceBuffer.col(replica).asVector() * deltaT;
        }
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    bool RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::checkCentroid() const {
        constexpr bool isGood = true;
        PositionMatrix centroid = ringPolymer.makeCentroidPos();
        cell.toDirect(centroid);
        for (auto& elem : centroid)
            if (!(PosScalarType(0) <= elem && elem <= PosScalarType(0)))
                return !isGood;
        return isGood;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, PosScalarType, Dim, NumReplica, ForceMatrixAllocator>::checkParam() const {
        const ScalarType cycle = PlainScalar(2 * M_PI) / ringPolymer.calcOmegaW(temperatureT);
        bool isSmallEnough = timeStep < cycle / PlainScalar(4);
        if (!isSmallEnough)
            throw std::invalid_argument("[Error]: Time step is too large");
    }
}
