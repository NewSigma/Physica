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

namespace Physica::Core {
    template<class ScalarType, class PosScalarType, unsigned int Dim>
    RPMD<ScalarType, PosScalarType, Dim>::RPMD(MDCellType cell_,
                                               size_t numReplica,
                                               size_t numContract,
                                               ScalarType temperatureT_,
                                               ScalarType timeStep_)
            : cell(std::move(cell_))
            , fftContract(numContract, 1)
            , timeStep(std::move(timeStep_)) {
        assert(0 < numContract && numContract <= numReplica);
        ringPolymer = RingPolymer(cell, numReplica);

        const size_t dof = getDOF();
        forceBuffer.resize(dof, numReplica);
        if (isContractEnabled()) {
            posContract.resize(dof, numContract);
            forceContract.resize(dof, numContract);
        }

        setTemperature(temperatureT_);

        dynamicStep = DynamicStepImpl(ringPolymer.calcRepBeta(temperatureT), getKSpaceSize(), getNumReplica());
        checkParam();
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    RPMD<ScalarType, PosScalarType, Dim>&
    RPMD<ScalarType, PosScalarType, Dim>::operator=(RPMD<ScalarType, PosScalarType, Dim> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class ForceModel, class Executor>
    void RPMD<ScalarType, PosScalarType, Dim>::updateForce(const ForceModel& model) {
        if (isContractEnabled()) {
            auto kernel_short = [&](unsigned int replica) {
                MDCellType cell = phaseToCell(replica);
                cell.normalize();
                auto saveTo = forceBuffer.col(replica);
                saveTo = model.template force_short<Executor>(std::move(cell));
            };
            auto future_short = Executor::parallel_for(kernel_short, getNumReplica(), Executor::getNumThread());

            contract();
            auto kernel_long = [&](unsigned int contract) {
                MDCellType cell = contractToCell(contract);
                cell.normalize();
                auto saveTo = forceContract.col(contract);
                saveTo = model.template force_long<Executor>(std::move(cell));
            };
            auto future_long = Executor::parallel_for(kernel_long, getNumContract(), Executor::getNumThread());

            Executor::auto_wait(future_long);
            decontract();
            Executor::auto_wait(future_short);
        }
        else {
            auto kernel = [&](unsigned int replica) {
                MDCellType cell = phaseToCell(replica);
                cell.normalize();
                auto saveTo = forceBuffer.col(replica);
                saveTo = model.template force<Executor>(std::move(cell));
            };
            auto future = Executor::parallel_for(kernel, getNumReplica(), Executor::getNumThread());
            Executor::auto_wait(future);
        }
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class ForceModel, class Executor>
    void RPMD<ScalarType, PosScalarType, Dim>::nve_step(const ForceModel& model) {
        forceStep(timeStep * 0.5);
        dynamicStep.nve_step(ringPolymer, timeStep);
        updateForce<ForceModel, Executor>(model);
        forceStep(timeStep * 0.5);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class ForceModel, class Executor>
    void RPMD<ScalarType, PosScalarType, Dim>::nve_step_for(ScalarType duration, const ForceModel& force) {
        uint64_t step = double(duration / timeStep) + 0.5;
        for (uint64_t _ = 0; _ < step; ++_)
            nve_step<ForceModel, Executor>(force);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class Thermostat, class RandomGenerator, class ForceModel, class Executor>
    void RPMD<ScalarType, PosScalarType, Dim>::nvt_step(const Thermostat& thermostat, RandomGenerator& gen, const ForceModel& model) {
        forceStep(timeStep * 0.5);
        dynamicStep.nve_step(ringPolymer, timeStep * 0.5);
        thermostat.step(ringPolymer, gen, timeStep);
        dynamicStep.nve_step(ringPolymer, timeStep * 0.5);
        updateForce<ForceModel, Executor>(model);
        forceStep(timeStep * 0.5);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class Thermostat, class RandomGenerator, class ForceModel, class Executor>
    void RPMD<ScalarType, PosScalarType, Dim>::nvt_step_for(ScalarType duration, const Thermostat& thermostat, RandomGenerator& gen, const ForceModel& model) {
        uint64_t step = double(duration / timeStep) + 0.5;
        for (uint64_t _ = 0; _ < step; ++_)
            nvt_step<Thermostat, RandomGenerator, ForceModel, Executor>(thermostat, gen, model);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class Thermostat, class RandomGenerator, class Barostat, class ForceModel, class Executor>
    void RPMD<ScalarType, PosScalarType, Dim>::npt_step(const Thermostat& thermostat, RandomGenerator& gen, Barostat& barostat, const ForceModel& model) {
        barostat.forceStep(*this, timeStep * 0.5);
        dynamicStep.npt_Step(ringPolymer, cell, barostat, timeStep * 0.5);
        thermostat.step(ringPolymer, gen, timeStep);
        dynamicStep.npt_Step(ringPolymer, cell, barostat, timeStep * 0.5);
        updateForce<ForceModel, Executor>(model);
        barostat.forceStep(*this, timeStep * 0.5);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class Thermostat, class RandomGenerator, class Barostat, class ForceModel, class Executor>
    void RPMD<ScalarType, PosScalarType, Dim>::npt_step_for(
            ScalarType duration,
            const Thermostat& thermostat,
            RandomGenerator& gen,
            Barostat& barostat,
            const ForceModel& model) {
        uint64_t step = double(duration / timeStep) + 0.5;
        for (uint64_t _ = 0; _ < step; ++_)
            npt_step<Thermostat, RandomGenerator, Barostat, ForceModel, Executor>(thermostat, gen, barostat, model);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class RandomGenerator>
    void RPMD<ScalarType, PosScalarType, Dim>::initMomentum(RandomGenerator& gen) {
        return ringPolymer.initMomentum(temperatureT, gen);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::scaleVelocity() {
        ringPolymer.scaleVelocity(temperatureT);
    }
    /**
     * Carrying out this function every several steps may stable the simulation.
     */
    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::normalizeCentroid() {
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

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    typename RPMD<ScalarType, PosScalarType, Dim>::MDCellType RPMD<ScalarType, PosScalarType, Dim>::phaseToCell(size_t replica) const {
        return MDCellType(cell.getLattice(), ringPolymer.makeBeadPos(replica), cell.getMassVec());
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    typename RPMD<ScalarType, PosScalarType, Dim>::MDCellType RPMD<ScalarType, PosScalarType, Dim>::contractToCell(size_t contract) const {
        PositionMatrix pos(getNumParticle(), Dim);
        auto phase = posContract.col(contract);
        size_t index = 0;
        for (auto& elem : pos) {
            elem = PosScalarType(phase[index]);
            ++index;
        }
        return MDCellType(cell.getLattice(), std::move(pos), cell.getMassVec());
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    typename RPMD<ScalarType, PosScalarType, Dim>::MDCellType RPMD<ScalarType, PosScalarType, Dim>::makeAverageCell() const {
        return MDCellType(getLattice(), ringPolymer.makeCentroidPos(), cell.getMassVec());
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::checkParam() const {
        const ScalarType cycle = ScalarType(2 * M_PI) / dynamicStep.getOmegaW();
        bool isSmallEnough = timeStep < cycle / ScalarType(4);
        if (!isSmallEnough)
            throw std::invalid_argument("[Error]: Time step is too large");
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::swap(RPMD& obj) noexcept {
        cell.swap(obj.cell);
        ringPolymer.swap(obj.ringPolymer);
        forceBuffer.swap(obj.forceBuffer);

        fftContract.swap(obj.fftContract);
        posContract.swap(obj.posContract);
        forceContract.swap(obj.forceContract);

        dynamicStep.swap(obj.dynamicStep);

        temperatureT.swap(obj.temperatureT);
        timeStep.swap(obj.timeStep);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class ForceModel>
    ScalarType RPMD<ScalarType, PosScalarType, Dim>::getClassicalPotentialEnergy(const ForceModel& model) const {
        ScalarType result = 0;
        for (size_t i = 0; i < getNumReplica(); ++i)
            result += model.potentialEnergy(phaseToCell(i));
        return result;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    ScalarType RPMD<ScalarType, PosScalarType, Dim>::getClassicalElastic() const {
        const size_t dof = getDOF();
        auto pos = getPhaseMatrix().bottomRows(dof);
        ScalarType result = 0;
        for (size_t i = 0; i < dof; ++i) {
            const ScalarType mass = cell.getMass(i / Dim);
            for (size_t j = 0; j < getNumReplica(); ++j)
                result += mass * square(dynamicStep.getOmegaW() * (pos(i, j) - pos(i, (j + 1) % getNumReplica()))) * 0.5;
        }
        return result;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class ForceModel>
    ScalarType RPMD<ScalarType, PosScalarType, Dim>::getClassicalInternalEnergy(const ForceModel& model) const {
        return ringPolymer.getClassicalKinetic() + getClassicalPotentialEnergy<ForceModel>(model) + getClassicalElastic();
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    ScalarType RPMD<ScalarType, PosScalarType, Dim>::calcKinetic() const {
        const ScalarType repBeta = ringPolymer.calcRepBeta(temperatureT);
        const size_t dof = getDOF();
        Vector<ScalarType> averaged_pos(dof, 0);
        for (size_t i = 0; i < dof; ++i)
            averaged_pos[i] = mean(getPhaseMatrix().row(dof + i));

        ScalarType kinetic = repBeta * dof;
        for (size_t replica = 0; replica < getNumReplica(); ++replica) {
            auto phase = getPhaseMatrix().col(replica);
            auto pos = phase.tail(dof);
            kinetic += (averaged_pos - pos) * forceBuffer.col(replica);
        }
        kinetic /= ScalarType(getNumReplica() * 2);
        return kinetic;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class ForceModel, class Executor>
    ScalarType RPMD<ScalarType, PosScalarType, Dim>::calcPotential(const ForceModel& model) const {
        Vector<ScalarType> temp(getNumReplica());
        auto kernel = [this, model, &temp](unsigned int replica) {
            temp[replica] = model.potentialEnergy(phaseToCell(replica));
        };
        Executor::parallel_for(kernel, getNumReplica(), Executor::getNumThread()).wait();
        return mean(temp);
    }
    /**
     * The function has size effect, extend the cell shall ease the problem, refer to [1].
     * 
     * Reference:
     * [1] M. J. Louwerse and E. J. Baerends, Chem. Phys. Lett. 421, 138 (2006).
     * [2] Thompson, Plimpton, Mattson, J Chem Phys, 131, 154107 (2009).
     */
    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class ForceModel>
    typename RPMD<ScalarType, PosScalarType, Dim>::LatticeMatrix RPMD<ScalarType, PosScalarType, Dim>::makeStress(const ForceModel& model) const {
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

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::setTemperature(ScalarType temperature) {
        assert(!temperature.isNegative());
        temperatureT = temperature;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::toContractBeadRepr(size_t posID) {
        assert(posID < getDOF());
        auto row = ringPolymer.getBuffer().row(1);
        auto head = row.head(fftContract.getKSpaceSize());
        fftContract.invTransform(head);
        auto pos = posContract.row(posID);
        pos = fftContract.getRSpace() * (ScalarType(getNumContract()) / ScalarType(getNumReplica()));
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::forceToNormRepr(size_t posID) {
        assert(posID < getDOF());
        fftContract.transform(forceContract.row(posID));
        auto row = ringPolymer.getBuffer().row(0);
        row = ScalarType(0);
        auto head = row.head(fftContract.getKSpaceSize());
        head = fftContract.getKSpace();
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::forceToBeadRepr(size_t posID) {
        assert(posID < getDOF());
        auto& buffer = ringPolymer.getBuffer();
        auto& fft = ringPolymer.getCanonicalFFT();
        fft.invTransform(buffer.row(0));
        auto f = forceBuffer.row(posID);
        f += fft.getRSpace() * (ScalarType(getNumReplica()) / ScalarType(getNumContract()));
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::contract() {
        const size_t dof = getDOF();
        for (size_t i = 0; i < dof; ++i) {
            ringPolymer.toNormalRepr(i);
            toContractBeadRepr(i);
        }
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::decontract() {
        const size_t dof = getDOF();
        for (size_t i = 0; i < dof; ++i) {
            forceToNormRepr(i);
            forceToBeadRepr(i);
        }
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::forceStep(ScalarType deltaT) {
        const size_t dof = getDOF();
        for (size_t replica = 0; replica < getNumReplica(); ++replica) {
            auto col = getPhaseMatrix().col(replica);
            auto momentum = col.head(dof);
            momentum += forceBuffer.col(replica).asVector() * deltaT;
        }
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    bool RPMD<ScalarType, PosScalarType, Dim>::checkCentroid() const {
        constexpr bool isGood = true;
        PositionMatrix centroid = ringPolymer.makeCentroidPos();
        cell.toDirect(centroid);
        for (auto& elem : centroid)
            if (!(PosScalarType::Zero() <= elem && elem <= PosScalarType::One()))
                return !isGood;
        return isGood;
    }
}
