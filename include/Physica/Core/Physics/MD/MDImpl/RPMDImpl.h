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
    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::RPMD(
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

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::read(const H5Location& loc, const char* name) {
        const auto group = loc.openGroup(name);
        LatticeMatrix temp{};
        temp.read(group, "lattice");
        setLattice(std::move(temp));
        getPhaseMatrix().read(group, "phase");
        forceBuffer.read(group, "force");
    }
    
    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::write(H5Location& loc, const char* name) const {
        auto group = loc.createGroup(name);
        getLattice().write(group, "lattice");
        getPhaseMatrix().write(group, "phase");
        forceBuffer.write(group, "force");
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>&
    RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::operator=(
            RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator> obj) noexcept {
        swap(obj);
        return *this;
    }
    /**
     * Contract method to improve performance in [1]
     * Reference:
     * [1] T. E. Markland, D. E. Manolopoulos. An efficient ring polymer contraction scheme for imaginary time path integral simulations[J]. J. Chem. Phys. 129, 024105 (2008)
     */
    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class ForceModel, class Executor>
    void RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::updateForce(ForceModel& model) {
        constexpr bool IsPeriodBoundary = Internal::Traits<ForceModel>::IsPeriodBoundary;
        constexpr bool isCudaEnabled = Internal::Traits<Executor>::isCudaEnabled;
        static_assert(!isCudaEnabled || std::allocator_traits<ForceMatrixAllocator>::isPageLocked
                , "[Error]: Allocator is not page locked, performance will decrease");
        if (isContractEnabled()) {
            auto kernel_short = [&](unsigned int thread) {
                const auto range = Executor::splitJob(getNumReplica(), Executor::getNumThread(), thread);
                for (size_t replica = range.first; replica < range.second; ++replica) {
                    MDCellType cell = phaseToCell(replica);
                    if constexpr (IsPeriodBoundary)
                        cell.normalize();
                    auto saveTo = forceBuffer.col(replica);
                    saveTo = model.template force_short<Executor>(std::move(cell));
                }
            };
            auto future_short = Executor::parallel_for(kernel_short, Executor::getNumThread());

            contract();
            auto kernel_long = [&](unsigned int thread) {
                const auto range = Executor::splitJob(getNumContract(), Executor::getNumThread(), thread);
                for (size_t contract = range.first; contract < range.second; ++contract) {
                    MDCellType cell = contractToCell(contract);
                    if constexpr (IsPeriodBoundary)
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
                if constexpr (IsPeriodBoundary)
                    cell.normalize();
                auto saveTo = forceBuffer.col(replica);
                using VectorType = typename decltype(saveTo)::VectorBase;
                model.template forceAsync<VectorType, Executor>(std::move(cell), saveTo);
            };
            auto future = Executor::parallel_for(kernel, getNumReplica());
            Executor::auto_wait(future);
            Executor::wait();
        }
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel,
             class ForceModel,
             class Executor>
    void RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::nve_step(KineticModel& kineticModel, ForceModel& forceModel) {
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        constexpr bool IsPeriodBoundary1 = Internal::Traits<KineticModel>::IsPeriodBoundary;
        constexpr bool IsPeriodBoundary2 = Internal::Traits<ForceModel>::IsPeriodBoundary;
        static_assert(isFreeModel || (IsPeriodBoundary1 == IsPeriodBoundary2), "[Error]: Inconsistent boundary condition");
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

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel,
             class ForceModel,
             class Executor>
    void RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::nve_step_for(
            ScalarType duration, KineticModel& kineticModel, ForceModel& forceModel) {
        const uint64_t step = Base::durationToStep(duration, timeStep);
        for (uint64_t _ = 0; _ < step; ++_)
            nve_step<KineticModel, ForceModel, Executor>(kineticModel, forceModel);
    }
    /**
     * BAOAB integrator as introduced in [1]
     * 
     * Reference:
     * [1] J. Chem. Phys, 145, 024103 (2016); https://doi.org/10.1063/1.4954990
     */
    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class Thermostat,
             class RandomPoolType,
             class KineticModel,
             class ForceModel,
             class Executor>
    void RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::nvt_step(
            const Thermostat& thermostat,
            KineticModel& kineticModel,
            ForceModel& forceModel) {
        constexpr bool isSeedFixed = RandomPoolType::isSeedFixed();
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        constexpr bool IsPeriodBoundary1 = Internal::Traits<KineticModel>::IsPeriodBoundary;
        constexpr bool IsPeriodBoundary2 = Internal::Traits<ForceModel>::IsPeriodBoundary;
        static_assert(isFreeModel || (IsPeriodBoundary1 == IsPeriodBoundary2), "[Error]: Inconsistent boundary condition");
        using NoRandExecutor = typename std::conditional<isSeedFixed, SequentialExecutor, Executor>::type;
        if (isFreeModel) {
            kineticModel.nve_step(ringPolymer, timeStep * 0.5);
            thermostat.template step<RandomPoolType, NoRandExecutor>(ringPolymer, timeStep);
            kineticModel.nve_step(ringPolymer, timeStep * 0.5);
        }
        else {
            forceStep(timeStep * 0.5);
            kineticModel.nve_step(ringPolymer, timeStep * 0.5);
            thermostat.template step<RandomPoolType, NoRandExecutor>(ringPolymer, timeStep);
            kineticModel.nve_step(ringPolymer, timeStep * 0.5);
            updateForce<ForceModel, Executor>(forceModel);
            forceStep(timeStep * 0.5);
        }
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class Thermostat,
             class RandomPoolType,
             class KineticModel,
             class ForceModel,
             class Executor>
    void RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::nvt_step_for(
            ScalarType duration,
            const Thermostat& thermostat,
            KineticModel& kineticModel,
            ForceModel& forceModel) {
        const uint64_t step = Base::durationToStep(duration, timeStep);
        for (uint64_t _ = 0; _ < step; ++_)
            nvt_step<Thermostat, RandomPoolType, KineticModel, ForceModel, Executor>(thermostat, kineticModel, forceModel);
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class Thermostat,
             class RandomPoolType,
             class Barostat,
             class KineticModel,
             class ForceModel,
             class Executor>
    void RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::npt_step(
            const Thermostat& thermostat,
            Barostat& barostat,
            KineticModel& kineticModel,
            ForceModel& forceModel) {
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        constexpr bool IsPeriodBoundary1 = Internal::Traits<KineticModel>::IsPeriodBoundary;
        constexpr bool IsPeriodBoundary2 = Internal::Traits<ForceModel>::IsPeriodBoundary;
        static_assert(isFreeModel || (IsPeriodBoundary1 == IsPeriodBoundary2), "[Error]: Inconsistent boundary condition");

        constexpr unsigned int BarostatOrder = Internal::Traits<Barostat>::Order;
        static_assert(BarostatOrder == 2 || BarostatOrder == 1, "[Error]: Invalid barostat");
        constexpr bool isSeedFixed = RandomPoolType::isSeedFixed();
        using NoRandExecutor = typename std::conditional<isSeedFixed, SequentialExecutor, Executor>::type;
        if constexpr (BarostatOrder == 2) {
            barostat.forceStep(*this, forceModel, timeStep * 0.5);
            kineticModel.npt_step(*this, barostat, timeStep * 0.5);
            thermostat.template step<RandomPoolType, NoRandExecutor>(ringPolymer, timeStep);
            kineticModel.npt_step(*this, barostat, timeStep * 0.5);
            updateForce<ForceModel, Executor>(forceModel);
            barostat.forceStep(*this, forceModel, timeStep * 0.5);
        }
        else {
            const LatticeMatrix stress = makeStressPrim<ForceModel, Executor>(forceModel);
            forceStep(timeStep * 0.5);
            kineticModel.nve_step(ringPolymer, timeStep * 0.5);
            barostat.template npt_step<This, ForceModel>(*this, stress, timeStep * 0.5);
            thermostat.template step<RandomPoolType, NoRandExecutor>(ringPolymer, timeStep);
            barostat.template npt_step<This, ForceModel>(*this, stress, timeStep * 0.5);
            kineticModel.nve_step(ringPolymer, timeStep * 0.5);
            if constexpr (Internal::Traits<ForceModel>::IsLatticeDependent)
                forceModel.setLattice(getLattice());
            updateForce<ForceModel, Executor>(forceModel);
            forceStep(timeStep * 0.5);
        }
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class Thermostat, class RandomPoolType, class Barostat, class KineticModel, class ForceModel, class Executor>
    void RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::npt_step_for(
            ScalarType duration,
            const Thermostat& thermostat,
            Barostat& barostat,
            KineticModel& kineticModel,
            ForceModel& forceModel) {
        const uint64_t step = Base::durationToStep(duration, timeStep);
        for (uint64_t _ = 0; _ < step; ++_)
            npt_step<Thermostat, RandomPoolType, Barostat, KineticModel, ForceModel, Executor>(thermostat, barostat, kineticModel, forceModel);
    }
    /**
     * Euler Semi-implicit integrator as introduced in [1]
     * 
     * Reference:
     * [1] Comput. Mater. Sci. 175, 109584 (2020); https://doi.org/10.1016/j.commatsci.2020.109584
     */
    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel, class ForceModel, class Executor>
    void RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::fire_vstep(
            FireModelType& fire, KineticModel& kineticModel, ForceModel& forceModel) {
        constexpr bool IsPeriodBoundary1 = Internal::Traits<KineticModel>::IsPeriodBoundary;
        constexpr bool IsPeriodBoundary2 = Internal::Traits<ForceModel>::IsPeriodBoundary;
        static_assert(IsPeriodBoundary1 == IsPeriodBoundary2, "[Error]: Inconsistent boundary condition");
        static_assert(NumReplica == 1, "[Error]: Relaxing using PIMD makes no sence, NumReplica = 1 shall be enough");
        static_assert(!Internal::is_empty_force_model<ForceModel>::value, "[Error]: Relax a empty model does nothing");

        kineticModel.nve_step(ringPolymer, timeStep);
        updateForce<ForceModel, Executor>(forceModel);
        fire.paramStep(*this);
        forceStep(fire.getTimeStep());
        fire.mixingStep(*this);
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<BaroType Type, class KineticModel, class ForceModel, class Executor>
    void RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::fire_pstep(
            FireModelType& fire,
            Berendsen<ScalarType, NumReplica, Type>& barostat,
            KineticModel& kineticModel,
            ForceModel& forceModel) {
        constexpr bool IsPeriodBoundary1 = Internal::Traits<KineticModel>::IsPeriodBoundary;
        constexpr bool IsPeriodBoundary2 = Internal::Traits<ForceModel>::IsPeriodBoundary;
        static_assert(IsPeriodBoundary1 == IsPeriodBoundary2, "[Error]: Inconsistent boundary condition");
        static_assert(NumReplica == 1, "[Error]: Relaxing using PIMD makes no sence, NumReplica = 1 shall be enough");
        static_assert(!Internal::is_empty_force_model<ForceModel>::value, "[Error]: Relax a empty model does nothing");

        kineticModel.nve_step(ringPolymer, timeStep);
        const LatticeMatrix stress = forceModel.virial(phaseToCell(0));
        barostat.template npt_step<ForceModel>(*this, stress, timeStep);
        if constexpr (Internal::Traits<ForceModel>::IsLatticeDependent)
            forceModel.setLattice(getLattice());
        updateForce<ForceModel, Executor>(forceModel);
        fire.paramStep(*this);
        forceStep(fire.getTimeStep());
        fire.mixingStep(*this);
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class RandomGenerator>
    void RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::initMomentum(RandomGenerator& gen) {
        return ringPolymer.initMomentum(temperatureT, gen);
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::scaleVelocity() {
        ringPolymer.scaleVelocity(temperatureT);
    }
    /**
     * Carrying out this function every several steps may stable the simulation.
     */
    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::normalizeCentroid() {
        PositionMatrix centroid = ringPolymer.makeCentroidPos();
        cell.toDirect(centroid);
        size_t index = getDOF();
        for (const auto elem : centroid) {
            const size_t component = index % Dim;
            const size_t atom_start = index - component;
            const int integer = float(elem);
            const Vector<ScalarType, Dim> delta = ScalarType(integer - elem.isNegative()) * cell.getLattice().row(component).asVector();
            for (size_t i = 0; i < Dim; ++i) {
                auto row = getPhaseMatrix().row(atom_start + i);
                row.asVector() -= delta[i];
            }
            ++index;
        }
        assert(checkCentroid());
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    typename RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::MDCellType
    RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::phaseToCell(size_t replica) const {
        return MDCellType(cell.getLattice(), ringPolymer.makeBeadPos(replica), cell.getMassVec());
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    typename RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::MDCellType
    RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::contractToCell(size_t contract) const {
        assert(contract < getNumContract());
        PositionMatrix pos(getNumParticle(), Dim);
        auto phase = posContract.col(contract);
        size_t index = 0;
        for (auto& elem : pos) {
            elem = ScalarType(phase[index]);
            ++index;
        }
        return MDCellType(cell.getLattice(), std::move(pos), cell.getMassVec());
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    typename RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::MDCellType
    RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::makeAverageCell() const {
        return MDCellType(getLattice(), ringPolymer.makeCentroidPos(), cell.getMassVec());
    }
    /**
     * Kinetic energy using virial estimator referenced from [1]
     * 
     * Reference:
     * [1] M. F. Herman, E. J. Bruskin, and B. J. Berne, J. Chem. Phys. 76, 5150(1982).
     */
    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    ScalarType RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::calcKinetic() const {
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

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    ScalarType RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::calcKinetic(size_t dofIndex) const {
        const ScalarType repBeta = ringPolymer.calcRepBeta(ringPolymer.calcTemperature());
        const auto pos = getPhaseMatrix().row(getDOF() + dofIndex);
        const ScalarType centroidPos = mean(pos);

        ScalarType kinetic = repBeta - (pos.asVector() - centroidPos) * forceBuffer.row(dofIndex).asVector();
        kinetic /= ScalarType(getNumReplica() * 2);
        return kinetic;
    }
    /**
     * Kinetic energy using primitive estimator referenced from [1]
     * Note: Use this estimator if NumReplica is small or if force model is \class EmptyForceModel
     * 
     * Reference:
     * [1] M. F. Herman, E. J. Bruskin, and B. J. Berne, J. Chem. Phys. 76, 5150 (1982); https://doi.org/10.1063/1.442815
     */
    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    ScalarType RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::calcKineticPrim() const {
        ScalarType kinetic = 0;
        for (size_t i = 0; i < getDOF(); ++i)
            kinetic += calcKineticPrim(i);
        return kinetic;
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    ScalarType RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::calcKineticPrim(size_t dofIndex) const {
        const ScalarType repBeta = ringPolymer.calcRepBeta(ringPolymer.calcTemperature());
        const ScalarType omegaW = ringPolymer.calcOmegaW(temperatureT);
        const auto pos = getPhaseMatrix().row(getDOF() + dofIndex);

        ScalarType kinetic = 0;
        for (size_t i = 0; i < getNumReplica(); ++i)
            kinetic += square(pos[i] - pos[(i + 1) % getNumReplica()]);
        const ScalarType factor = getMassVec()[dofIndex / Dim] * square(omegaW) / ScalarType(getNumReplica());
        kinetic = (-kinetic * factor + repBeta) * PlainScalar(0.5);
        return kinetic;
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    inline ScalarType RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::calcKineticClassical() const {
        return ringPolymer.calcKineticClassical() / square(ScalarType(getNumReplica()));
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class ForceModel, class Executor>
    ScalarType RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::calcPotential(const ForceModel& model) const {
        Vector<ScalarType> temp(getNumReplica());
        auto kernel = [this, model, &temp](unsigned int replica) {
            temp[replica] = model.potentialEnergy(phaseToCell(replica));
        };
        Executor::parallel_for(kernel, getNumReplica(), Executor::getNumThread()).wait();
        return mean(temp);
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class ForceModel>
    ScalarType RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::calcPotentialClassical(const ForceModel& model) const {
        ScalarType result = 0;
        for (size_t i = 0; i < getNumReplica(); ++i)
            result += model.potentialEnergy(phaseToCell(i));
        return result;
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    ScalarType RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::calcClassicalElastic() const {
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

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class ForceModel>
    ScalarType RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::calcClassicalInternalEnergy(const ForceModel& model) const {
        return calcKineticClassical() + calcPotentialClassical<ForceModel>(model) + calcClassicalElastic();
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class ForceModel, class Executor>
    ScalarType RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::calcPressThermo(ForceModel& model) const {
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        constexpr bool IsPeriodBoundary = Internal::Traits<ForceModel>::IsPeriodBoundary;
        ScalarType result = calcKinetic() / (getVolume() * (Dim * 0.5));
        if constexpr (!isFreeModel) {
            const size_t numReplica = getNumReplica();
            ScalarType temp = 0;
            for (size_t i = 0; i < numReplica; ++i) {
                auto cell = phaseToCell(i);
                if constexpr (IsPeriodBoundary)
                    cell.normalize();
                temp += model.virial(std::move(cell)).trace();
            }
            result += temp / ScalarType(numReplica * Dim);
        }
        return result;
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class ForceModel, class Executor>
    typename RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::LatticeMatrix
    RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::makeStressPrim(ForceModel& model) const {
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        if constexpr (NumReplica == 1)
            return makeStressClassical<ForceModel, Executor>(model);

        Utils::Array<LatticeMatrix> buffer(getNumReplica());
        const ScalarType squaredOmegaW = square(ringPolymer.calcOmegaW(temperatureT));
        auto kernel = [this, &model, &buffer, squaredOmegaW](unsigned int replica) {
            using VectorType = Vector<ScalarType, Dim>;
            const size_t dof = getDOF();
            const size_t numReplica = getNumReplica();
            const auto col = getPhaseMatrix().col(replica);
            const auto momentum = col.head(dof);
            const auto pos = col.tail(dof);
            const auto col1 = getPhaseMatrix().col((replica + 1) % numReplica);
            const auto pos1 = col1.tail(dof);

            LatticeMatrix kineticStress(Dim, Dim, 0);
            for (size_t i = 0; i < getNumParticle(); ++i) {
                const size_t from = i * Dim;
                const size_t to = from + Dim;
                const ScalarType mass = getMassVec()[i];
                const ScalarType repMass = reciprocal(mass);
                const auto atomMomentum = momentum.segment(from, to);
                kineticStress += (repMass * atomMomentum) * atomMomentum.transpose();

                const auto atomPos = pos.segment(from, to);
                const auto atomPos1 = pos1.segment(from, to);
                const VectorType deltaPos = atomPos - atomPos1;
                const ScalarType factorK = mass * squaredOmegaW;
                kineticStress -= (factorK * deltaPos) * deltaPos.transpose();
            }
            buffer[replica] = kineticStress * reciprocal(getVolume());

            if constexpr (!isFreeModel) {
                constexpr bool IsPeriodBoundary = Internal::Traits<ForceModel>::IsPeriodBoundary;
                MDCellType cell = phaseToCell(replica);
                if constexpr (IsPeriodBoundary)
                    cell.normalize();
                buffer[replica] += model.virial(cell);
            }
        };
        Executor::parallel_for(kernel, getNumReplica(), Executor::getNumThread()).wait();
        LatticeMatrix result(Dim, Dim, 0);
        for (size_t i = 0; i < buffer.getLength(); ++i)
            toNextMean(result, i, buffer[i]);
        return result;
    }
    /**
     * Reference:
     * [1] Comp. Phys. Comm. 185, 1019 (2013); https://doi.org/10.1016/j.cpc.2013.10.027
     */
    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class ForceModel, class Executor>
    typename RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::LatticeMatrix
    RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::makeStressCentroid(ForceModel& model) const {
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        static_assert(!isFreeModel, "[Error]: This function does not apply to ideal gas model");
        if constexpr (NumReplica == 1)
            return makeStressClassical<ForceModel, Executor>(model);

        Utils::Array<LatticeMatrix> buffer(getNumReplica());
        const auto centroidPos = ringPolymer.makeCentroidPos();
        auto kernel = [this, &model, &buffer, &centroidPos](unsigned int replica) {
            using VectorType = Vector<ScalarType, Dim>;
            const size_t dof = getDOF();
            const size_t numReplica = getNumReplica();
            const auto col = getPhaseMatrix().col(replica);
            const auto momentum = col.head(dof);
            const auto pos = col.tail(dof);
            const auto centroid = centroidPos.flatten();
            const auto force = forceBuffer.col(replica);

            LatticeMatrix classicalKineticStress(Dim, Dim, 0);
            LatticeMatrix quantumKineticStress(Dim, Dim, 0);
            for (size_t i = 0; i < getNumParticle(); ++i) {
                const size_t from = i * Dim;
                const size_t to = from + Dim;
                const ScalarType mass = getMassVec()[i];
                const ScalarType factor = reciprocal(mass * ScalarType(numReplica));
                const auto atomMomentum = momentum.segment(from, to);
                classicalKineticStress += (factor * atomMomentum) * atomMomentum.transpose();

                const auto atomPos = pos.segment(from, to);
                const auto atomCentroid = centroid.segment(from, to);
                const auto atomForce = force.segment(from, to);
                const VectorType deltaPos = atomPos - atomCentroid;
                quantumKineticStress += deltaPos * atomForce.transpose();
            }
            buffer[replica] = (classicalKineticStress - quantumKineticStress * ScalarType(0.5)) * reciprocal(getVolume());

            if constexpr (!isFreeModel) {
                constexpr bool IsPeriodBoundary = Internal::Traits<ForceModel>::IsPeriodBoundary;
                MDCellType cell = phaseToCell(replica);
                if constexpr (IsPeriodBoundary)
                    cell.normalize();
                buffer[replica] += model.virial(cell);
            }
        };
        Executor::parallel_for(kernel, getNumReplica(), Executor::getNumThread()).wait();
        LatticeMatrix result(Dim, Dim, 0);
        for (size_t i = 0; i < buffer.getLength(); ++i)
            toNextMean(result, i, buffer[i]);
        return result;
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class ForceModel, class Executor>
    typename RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::LatticeMatrix
    RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::makeStressVirial(ForceModel& model) const {
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        static_assert(!isFreeModel, "[Error]: This function does not apply to ideal gas model");
        if constexpr (NumReplica == 1)
            return makeStressClassical<ForceModel, Executor>(model);

        Utils::Array<LatticeMatrix> buffer(getNumReplica());
        const auto centroidPos = ringPolymer.makeCentroidPos();
        auto kernel = [this, &model, &buffer, &centroidPos](unsigned int replica) {
            using VectorType = Vector<ScalarType, Dim>;
            constexpr bool IsPeriodBoundary = Internal::Traits<ForceModel>::IsPeriodBoundary;
            const size_t dof = getDOF();
            const size_t numReplica = getNumReplica();
            const auto col = getPhaseMatrix().col(replica);
            const auto momentum = col.head(dof);
            const auto pos = col.tail(dof);
            const auto centroid = centroidPos.flatten();
            const auto force = forceBuffer.col(replica);
            auto cell = phaseToCell(replica);
            if constexpr (IsPeriodBoundary)
                cell.normalize();
            const auto forceConst = model.forceConst(cell);

            LatticeMatrix classicalKineticStress(Dim, Dim, 0);
            LatticeMatrix quantumKineticStress(Dim, Dim, 0);
            for (size_t i = 0; i < getNumParticle(); ++i) {
                const size_t from = i * Dim;
                const size_t to = from + Dim;
                const ScalarType mass = getMassVec()[i];
                const ScalarType factor = reciprocal(mass * ScalarType(numReplica));
                const auto atomMomentum = momentum.segment(from, to);
                classicalKineticStress += (factor * atomMomentum) * atomMomentum.transpose();

                const auto atomPos = pos.segment(from, to);
                const auto atomCentroid = centroid.segment(from, to);
                const auto atomForce = force.segment(from, to);
                const VectorType deltaPos = atomPos - atomCentroid;
                quantumKineticStress += deltaPos * atomForce.transpose();

                for (size_t j = 0; j < getNumParticle(); ++j) {
                    const auto atomPos1 = pos.segment(j * Dim, (j + 1) * Dim);
                    const auto block = forceConst.block(Dim * j, Dim, Dim * i, Dim);
                    const VectorType temp = block * deltaPos;
                    quantumKineticStress -= atomPos1 * temp.transpose();
                }
            }
            buffer[replica] = (classicalKineticStress - quantumKineticStress * ScalarType(0.5)) * reciprocal(getVolume());

            if constexpr (!isFreeModel)
                buffer[replica] += model.virial(cell);
        };
        Executor::parallel_for(kernel, getNumReplica(), Executor::getNumThread()).wait();
        LatticeMatrix result(Dim, Dim, 0);
        for (size_t i = 0; i < buffer.getLength(); ++i)
            toNextMean(result, i, buffer[i]);
        return result;
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class ForceModel, class Executor>
    typename RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::LatticeMatrix
    RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::makeStressClassical(ForceModel& model) const {
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        constexpr bool IsPeriodBoundary = Internal::Traits<ForceModel>::IsPeriodBoundary;
        LatticeMatrix result(Dim, Dim, 0);
        if constexpr (NumReplica == 1) {
            const ScalarType repVolume = reciprocal(getVolume());
            const auto col = getPhaseMatrix().col(0);
            const auto momentum = col.head(getDOF());
            for (size_t i = 0; i < getNumParticle(); ++i) {
                const size_t from = i * Dim;
                const size_t to = from + Dim;
                const ScalarType repMass = reciprocal(getMassVec()[i]);
                const auto atomMomentum = momentum.segment(from, to);
                result += (repMass * atomMomentum) * atomMomentum.transpose();
            }
            result *= repVolume;

            if constexpr (!isFreeModel) {
                MDCellType cell = phaseToCell(0);
                if constexpr (IsPeriodBoundary)
                    cell.normalize();
                result += model.virial(cell);
            }
        }
        else {
            Utils::Array<LatticeMatrix> buffer(getNumReplica());
            auto kernel = [this, &model, &buffer](unsigned int replica) {
                const size_t dof = getDOF();
                const size_t numReplica = getNumReplica();
                const auto col = getPhaseMatrix().col(replica);
                const auto momentum = col.head(dof);
                const auto pos = col.tail(dof);

                LatticeMatrix kineticStress(Dim, Dim);
                for (size_t i = 0; i < getNumParticle(); ++i) {
                    const size_t from = i * Dim;
                    const size_t to = from + Dim;
                    const ScalarType mass = getMassVec()[i];
                    const ScalarType factor = reciprocal(mass * ScalarType(numReplica));
                    const auto atomMomentum = momentum.segment(from, to);
                    kineticStress += (factor * atomMomentum) * atomMomentum.transpose();
                }
                buffer[replica] = kineticStress * reciprocal(getVolume());

                if constexpr (!isFreeModel) {
                    MDCellType cell = phaseToCell(replica);
                    if constexpr (IsPeriodBoundary)
                        cell.normalize();
                    buffer[replica] += model.virial(cell);
                }
            };
            Executor::parallel_for(kernel, getNumReplica(), Executor::getNumThread()).wait();
            for (size_t i = 0; i < buffer.getLength(); ++i)
                toNextMean(result, i, buffer[i]);
        }
        return result;
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::swap(RPMD& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        cell.swap(obj.cell);
        ringPolymer.swap(obj.ringPolymer);
        forceBuffer.swap(obj.forceBuffer);

        fftContract.swap(obj.fftContract);
        posContract.swap(obj.posContract);
        forceContract.swap(obj.forceContract);

        temperatureT.swap(obj.temperatureT);
        timeStep.swap(obj.timeStep);
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::setTemperature(ScalarType temperature) {
        assert(!temperature.isNegative() && "[Error]: Negative temperature is not physical");
        temperatureT = temperature;
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::toContractBeadRepr(size_t posID) {
        assert(posID < getDOF());
        auto row = ringPolymer.getBuffer().row(1);
        auto head = row.head(fftContract.getKSpaceSize());
        fftContract.invTransform(head);
        auto pos = posContract.row(posID);
        pos = fftContract.getRSpace() * (ScalarType(getNumContract()) / ScalarType(getNumReplica()));
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::forceToNormRepr(size_t posID) {
        assert(posID < getDOF());
        fftContract.transform(forceContract.row(posID));
        auto row = ringPolymer.getBuffer().row(0);
        row.asVector() = ScalarType(0);
        auto head = row.head(fftContract.getKSpaceSize());
        head = fftContract.getKSpace();
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::forceToBeadRepr(size_t posID) {
        assert(posID < getDOF());
        auto& buffer = ringPolymer.getBuffer();
        auto& fft = ringPolymer.getCanonicalFFT();
        fft.invTransform(buffer.row(0));
        auto f = forceBuffer.row(posID);
        f.asVector() += fft.getRSpace() * (ScalarType(getNumReplica()) / ScalarType(getNumContract()));
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::contract() {
        const size_t dof = getDOF();
        for (size_t i = 0; i < dof; ++i) {
            ringPolymer.toNormalRepr(i);
            toContractBeadRepr(i);
        }
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::decontract() {
        const size_t dof = getDOF();
        for (size_t i = 0; i < dof; ++i) {
            forceToNormRepr(i);
            forceToBeadRepr(i);
        }
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::forceStep(ScalarType deltaT) {
        const size_t dof = getDOF();
        for (size_t replica = 0; replica < getNumReplica(); ++replica) {
            auto col = getPhaseMatrix().col(replica);
            auto momentum = col.head(dof);
            momentum += forceBuffer.col(replica).asVector() * deltaT;
        }
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    bool RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::checkCentroid() const {
        constexpr bool isGood = true;
        PositionMatrix centroid = ringPolymer.makeCentroidPos();
        cell.toDirect(centroid);
        for (auto& elem : centroid)
            if (!(ScalarType(0) <= elem && elem <= ScalarType(0)))
                return !isGood;
        return isGood;
    }

    template<class ScalarType, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>::checkParam() const {
        if constexpr (NumReplica != 1) {
            if (getNumReplica() == 1)
                throw std::invalid_argument("[Warning]: Set tparam NumReplica = 1 may gain better performance");
            const ScalarType cycle = PlainScalar(2 * M_PI) / ringPolymer.calcOmegaW(temperatureT);
            bool isSmallEnough = timeStep < cycle / PlainScalar(4);
            if (!isSmallEnough)
                throw std::invalid_argument("[Error]: Time step is too large");
        }
    }
}
