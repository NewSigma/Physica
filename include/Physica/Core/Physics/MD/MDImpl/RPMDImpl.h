/*
 * Copyright 2022-2025 Weibo He.
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
#include "Physica/Core/Math/Random/Random.h"
#include "../RPMD.h"

namespace Physica {
    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::RPMD(
            MDCellType cell_,
            size_t numReplica,
            size_t numContract,
            T temperatureT_,
            T timeStep_)
            : cell(std::move(cell_))
            , fftContract(numContract, PlanFlag::Estimate)
            , timeStep(std::move(timeStep_)) {
        assert(0 < numContract && numContract <= numReplica);
        assert(NumReplica == Dynamic || NumReplica == numReplica);
        ringPolymer = RingPolymerType(cell, numReplica);

        const size_t dof = getDOF();
        forceBuffer.resize(dof, numReplica, 0);
        if (isContractEnabled()) {
            posContract.resize(dof, numContract);
            forceContract.resize(dof, numContract);
        }

        setTemperature(temperatureT_);
        checkParam();
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    RPMD<T, Dim, NumReplica, ForceMatrixAllocator>&
    RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::operator=(
            RPMD<T, Dim, NumReplica, ForceMatrixAllocator> obj) noexcept {
        swap(obj);
        return *this;
    }
    /**
     * Contract method to improve performance introduced in [1]
     * 
     * Reference:
     * [1] J. Chem. Phys. 129, 024105 (2008); https://doi.org/10.1063/1.2953308
     */
    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class ForceModel, class Executor>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::updateForce(ForceModel& model) {
        constexpr bool IsPeriodBoundary = Traits<ForceModel>::IsPeriodBoundary;
        constexpr bool UseCUDA = Traits<Executor>::UseCUDA;
        static_assert(!UseCUDA || std::allocator_traits<ForceMatrixAllocator>::isPageLocked
                , "[Error]: Allocator is not page locked, performance will decrease");
        if (!isContractEnabled()) {
            auto kernel = [this, &model](unsigned int replica) {
                MDCellType cell = phaseToCell(replica);
                if constexpr (IsPeriodBoundary)
                    cell.normalize();
                auto saveTo = forceBuffer.col(replica);
                model.template forceAsync<decltype(saveTo), Executor>(std::move(cell), saveTo);
            };

            if constexpr (NumReplica == 1)
                kernel(0);
            else {
                auto future = Executor::parallel_for(kernel, getNumReplica());
                Executor::auto_wait(future);
            }
            Executor::wait();
            return;
        }

        constexpr bool IsContractable = Traits<ForceModel>::IsContractable;
        if constexpr (IsContractable) {
            auto kernel_uncontract = [&](unsigned int thread) {
                const auto range = Executor::splitJob(getNumReplica(), Executor::getNumThread(), thread);
                for (size_t replica = range.first; replica < range.second; ++replica) {
                    MDCellType cell = phaseToCell(replica);
                    if constexpr (IsPeriodBoundary)
                        cell.normalize();
                    auto saveTo = forceBuffer.col(replica);
                    saveTo = model.template force_uncontract<Executor>(std::move(cell));
                }
            };
            auto future_uncontract = Executor::parallel_for(kernel_uncontract, Executor::getNumThread());

            contract();
            auto kernel_contract = [&](unsigned int thread) {
                const auto range = Executor::splitJob(getNumContract(), Executor::getNumThread(), thread);
                for (size_t contract = range.first; contract < range.second; ++contract) {
                    MDCellType cell = contractToCell(contract);
                    if constexpr (IsPeriodBoundary)
                        cell.normalize();
                    auto saveTo = forceContract.col(contract);
                    saveTo = model.template force_contract<Executor>(std::move(cell));
                }
            };
            auto future_contract = Executor::parallel_for(kernel_contract, Executor::getNumThread());

            Executor::auto_wait(future_uncontract);
            Executor::auto_wait(future_contract);
            decontract();
        }
        else
            noImpl("[Error]: Force contract is not implemented");
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel,
             class ForceModel,
             class Executor>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::nve_step(KineticModel& kineticModel, ForceModel& forceModel) {
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        constexpr bool IsPeriodBoundary1 = Traits<KineticModel>::IsPeriodBoundary;
        constexpr bool IsPeriodBoundary2 = Traits<ForceModel>::IsPeriodBoundary;
        static_assert(isFreeModel || (IsPeriodBoundary1 == IsPeriodBoundary2), "[Error]: Inconsistent boundary condition");
        if (isFreeModel) {
            kineticModel.nve_step(ringPolymer, timeStep);
        }
        else {
            forceStep(timeStep * Tv(0.5));
            kineticModel.nve_step(ringPolymer, timeStep);
            updateForce<ForceModel, Executor>(forceModel);
            forceStep(timeStep * Tv(0.5));
        }
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel,
             class ForceModel,
             class Executor>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::nve_step_for(
            T duration, KineticModel& kineticModel, ForceModel& forceModel) {
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
    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class Thermostat,
             RNG R,
             class KineticModel,
             class ForceModel,
             class Executor>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::nvt_step(
            const Thermostat& thermostat,
            KineticModel& kineticModel,
            ForceModel& forceModel) {
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        constexpr bool IsPeriodBoundary1 = Traits<KineticModel>::IsPeriodBoundary;
        constexpr bool IsPeriodBoundary2 = Traits<ForceModel>::IsPeriodBoundary;
        static_assert(isFreeModel || (IsPeriodBoundary1 == IsPeriodBoundary2), "[Error]: Inconsistent boundary condition");

        constexpr bool isSeedFixed = Traits<R>::IsSeedFixed;
        using NoRandExecutor = std::conditional<isSeedFixed, SeqExecutor, Executor>::type;

        constexpr bool IsCentroidCoupled = Traits<Thermostat>::IsCentroidCoupled;
        if constexpr (IsCentroidCoupled && IsPeriodBoundary1)
            assert(thermostat.isRemoveDriftEnabled() && "[Error]: Because the KineticModel has period boundary, thermostat should remove its effect on centroid");

        if constexpr (isFreeModel) {
            kineticModel.nve_step(ringPolymer, timeStep * 0.5);
            thermostat.template step<R, NoRandExecutor>(ringPolymer, timeStep);
            kineticModel.nve_step(ringPolymer, timeStep * 0.5);
        }
        else {
            forceStep(timeStep * 0.5);
            kineticModel.nve_step(ringPolymer, timeStep * 0.5);
            thermostat.template step<R, NoRandExecutor>(ringPolymer, timeStep);
            kineticModel.nve_step(ringPolymer, timeStep * 0.5);
            updateForce<ForceModel, Executor>(forceModel);
            forceStep(timeStep * 0.5);
        }
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class Thermostat,
             RNG R,
             class KineticModel,
             class ForceModel,
             class Executor>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::nvt_step_for(
            T duration,
            const Thermostat& thermostat,
            KineticModel& kineticModel,
            ForceModel& forceModel) {
        const uint64_t step = Base::durationToStep(duration, timeStep);
        for (uint64_t _ = 0; _ < step; ++_)
            nvt_step<Thermostat, R, KineticModel, ForceModel, Executor>(thermostat, kineticModel, forceModel);
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class Thermostat,
             RNG R,
             class Barostat,
             class KineticModel,
             class ForceModel,
             class Executor>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::npt_step(
            const Thermostat& thermostat,
            Barostat& barostat,
            KineticModel& kineticModel,
            ForceModel& forceModel) {
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        constexpr bool IsPeriodBoundary1 = Traits<KineticModel>::IsPeriodBoundary;
        constexpr bool IsPeriodBoundary2 = Traits<ForceModel>::IsPeriodBoundary;
        static_assert(isFreeModel || (IsPeriodBoundary1 == IsPeriodBoundary2), "[Error]: Inconsistent boundary condition");

        constexpr unsigned int BarostatOrder = Traits<Barostat>::Order;
        static_assert(BarostatOrder == 2 || BarostatOrder == 1, "[Error]: Invalid barostat");

        constexpr bool isSeedFixed = Traits<R>::IsSeedFixed;
        using NoRandExecutor = std::conditional<isSeedFixed, SeqExecutor, Executor>::type;

        constexpr bool IsCentroidCoupled = Traits<Thermostat>::IsCentroidCoupled;
        if constexpr (IsCentroidCoupled && IsPeriodBoundary1)
            assert(thermostat.isRemoveDriftEnabled() && "[Error]: Because the KineticModel has period boundary, thermostat should remove its effect on centroid");

        if constexpr (BarostatOrder == 2) {
            barostat.forceStep(*this, forceModel, timeStep * 0.5);
            kineticModel.npt_step(*this, barostat, timeStep * 0.5);
            thermostat.template step<R, NoRandExecutor>(ringPolymer, timeStep);
            kineticModel.npt_step(*this, barostat, timeStep * 0.5);
            updateForce<ForceModel, Executor>(forceModel);
            barostat.forceStep(*this, forceModel, timeStep * 0.5);
        }
        else {
            const LatticeMatrix stress = makeStressPrim<ForceModel, Executor>(forceModel);
            forceStep(timeStep * 0.5);
            kineticModel.nve_step(ringPolymer, timeStep * 0.5);
            barostat.template npt_step<This, ForceModel>(*this, stress, timeStep * 0.5);
            thermostat.template step<R, NoRandExecutor>(ringPolymer, timeStep);
            barostat.template npt_step<This, ForceModel>(*this, stress, timeStep * 0.5);
            kineticModel.nve_step(ringPolymer, timeStep * 0.5);
            if constexpr (Traits<ForceModel>::IsLatticeDependent)
                forceModel.setLattice(getLattice());
            updateForce<ForceModel, Executor>(forceModel);
            forceStep(timeStep * 0.5);
        }
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class Thermostat, RNG R, class Barostat, class KineticModel, class ForceModel, class Executor>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::npt_step_for(
            T duration,
            const Thermostat& thermostat,
            Barostat& barostat,
            KineticModel& kineticModel,
            ForceModel& forceModel) {
        const uint64_t step = Base::durationToStep(duration, timeStep);
        for (uint64_t _ = 0; _ < step; ++_)
            npt_step<Thermostat, R, Barostat, KineticModel, ForceModel, Executor>(thermostat, barostat, kineticModel, forceModel);
    }
    /**
     * fire_vstep is fire_v(olume)step
     * 
     * Using euler semi-implicit integrator as introduced in [1]
     * 
     * Reference:
     * [1] Comput. Mater. Sci. 175, 109584 (2020); https://doi.org/10.1016/j.commatsci.2020.109584
     */
    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel, class ForceModel, class Executor>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::fire_vstep(
            FireModelType& fire, KineticModel& kineticModel, ForceModel& forceModel) {
        constexpr bool IsPeriodBoundary1 = Traits<KineticModel>::IsPeriodBoundary;
        constexpr bool IsPeriodBoundary2 = Traits<ForceModel>::IsPeriodBoundary;
        static_assert(IsPeriodBoundary1 == IsPeriodBoundary2, "[Error]: Inconsistent boundary condition");
        static_assert(NumReplica == 1, "[Error]: Relaxing using PIMD makes no sence, NumReplica = 1 shall be enough");
        static_assert(!Internal::is_empty_force_model<ForceModel>::value, "[Error]: Relax a empty model does nothing");

        kineticModel.nve_step(ringPolymer, timeStep);
        updateForce<ForceModel, Executor>(forceModel);
        fire.paramStep(*this);
        forceStep(fire.getTimeStep());
        fire.mixingStep(*this);
    }
    /**
     * fire_pstep is fire_p(ress)step
     */
    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<BaroType Type, class KineticModel, class ForceModel, class Executor>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::fire_pstep(
            CFireModel<T, Dim, Type>& cfire,
            KineticModel& kineticModel,
            ForceModel& forceModel) {
        constexpr bool IsPeriodBoundary1 = Traits<KineticModel>::IsPeriodBoundary;
        constexpr bool IsPeriodBoundary2 = Traits<ForceModel>::IsPeriodBoundary;
        static_assert(IsPeriodBoundary1 == IsPeriodBoundary2, "[Error]: Inconsistent boundary condition");
        static_assert(NumReplica == 1, "[Error]: Relaxing using PIMD makes no sence, NumReplica = 1 shall be enough");
        static_assert(!Internal::is_empty_force_model<ForceModel>::value, "[Error]: Relax a empty model does nothing");

        kineticModel.nve_step(ringPolymer, timeStep);
        cfire.nve_step(*this);
        if constexpr (Traits<ForceModel>::IsLatticeDependent)
            forceModel.setLattice(getLattice());
        updateForce<ForceModel, Executor>(forceModel);
        const LatticeMatrix stress = forceModel.virial(phaseToCell(0));
        cfire.paramStep(*this);
        forceStep(cfire.getTimeStep());
        cfire.forceStep(*this, std::move(stress));
        cfire.mixingStep(*this);
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel, RNG R>
    inline void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::initMomentum() {
        return ringPolymer.template initMomentum<KineticModel, R>(temperatureT);
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel>
    inline void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::scaleVelocity() {
        ringPolymer.template scaleVelocity<KineticModel>(temperatureT);
    }
    /**
     * Carrying out this function every several steps may stable the simulation.
     */
    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::normalizeCentroid() {
        PositionMatrix centroid = ringPolymer.makeCentroidPos();
        cell.toDirect(centroid);
        size_t index = getDOF();
        for (const auto elem : centroid) {
            const size_t component = index % Dim;
            const size_t atom_start = index - component;
            const int integer = float(elem);
            const DenseVector<T, Dim> delta = T(integer - elem.isNegative()) * cell.getLattice().row(component);
            for (size_t i = 0; i < Dim; ++i) {
                auto row = getPhaseMatrix().row(atom_start + i);
                row -= delta[i];
            }
            ++index;
        }
        assert(checkCentroid());
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::MDCellType
    RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::phaseToCell(size_t replica) const {
        return MDCellType(cell.getLattice(), ringPolymer.makeBeadPos(replica), cell.getMassVec());
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::MDCellType
    RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::contractToCell(size_t contract) const {
        assert(contract < getNumContract());
        PositionMatrix pos(getNumParticle(), Dim);
        auto phase = posContract.col(contract);
        size_t index = 0;
        for (auto& elem : pos.asArray()) {
            elem = T(phase[index]);
            ++index;
        }
        return MDCellType(cell.getLattice(), std::move(pos), cell.getMassVec());
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::MDCellType
    RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::makeAverageCell() const {
        return MDCellType(getLattice(), ringPolymer.makeCentroidPos(), cell.getMassVec());
    }
    /**
     * Kinetic energy using virial estimator referenced from [1]
     * 
     * Reference:
     * [1] M. F. Herman, E. J. Bruskin, and B. J. Berne, J. Chem. Phys. 76, 5150(1982).
     */
    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel>
    T RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::calcKinetic() const {
        const T repBeta = ringPolymer.calcRepBeta(calcTemperature<KineticModel>());
        const size_t dof = getDOF();
        const auto centroidPos = ringPolymer.makeCentroidPos();

        T kinetic = repBeta * dof;
        for (size_t replica = 0; replica < getNumReplica(); ++replica) {
            auto phase = getPhaseMatrix().col(replica);
            auto pos = phase.tail(dof);
            kinetic += (centroidPos.flatten() - pos) * forceBuffer.col(replica);
        }
        kinetic /= T(getNumReplica() * 2);
        return kinetic;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel>
    T RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::calcKinetic(size_t dofIndex) const {
        const T repBeta = ringPolymer.calcRepBeta(calcTemperature<KineticModel>());
        const auto pos = getPhaseMatrix().row(getDOF() + dofIndex);
        const T centroidPos = pos.mean();

        T kinetic = repBeta - (pos - centroidPos) * forceBuffer.row(dofIndex);
        kinetic /= T(getNumReplica() * 2);
        return kinetic;
    }
    /**
     * Kinetic energy using primitive estimator referenced from [1]
     * Note: Use this estimator if NumReplica is small or if force model is \class EmptyForceModel
     * 
     * Reference:
     * [1] M. F. Herman, E. J. Bruskin, and B. J. Berne, J. Chem. Phys. 76, 5150 (1982); https://doi.org/10.1063/1.442815
     */
    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel>
    T RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::calcKineticPrim() const {
        T kinetic = 0;
        for (size_t i = 0; i < getDOF(); ++i)
            kinetic += calcKineticPrim<KineticModel>(i);
        return kinetic;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel>
    T RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::calcKineticPrim(size_t dofIndex) const {
        const T repBeta = ringPolymer.calcRepBeta(calcTemperature<KineticModel>());
        const T omegaW = ringPolymer.calcOmegaW(temperatureT);
        const auto pos = getPhaseMatrix().row(getDOF() + dofIndex);

        T kinetic = 0;
        for (size_t i = 0; i < getNumReplica(); ++i)
            kinetic += square(pos[i] - pos[(i + 1) % getNumReplica()]);
        const T factor = getMassVec()[dofIndex / Dim] * square(omegaW) / T(getNumReplica());
        kinetic = (-kinetic * factor + repBeta) * Tv(0.5);
        return kinetic;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    inline T RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::calcKineticClassical() const {
        return ringPolymer.calcKineticClassical();
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class ForceModel, class Executor>
    T RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::calcPotential(const ForceModel& model) const {
        VectorND<T> temp(getNumReplica());
        auto kernel = [this, model, &temp](unsigned int replica) {
            temp[replica] = model.potentialV(phaseToCell(replica));
        };
        Executor::parallel_for(kernel, getNumReplica(), Executor::getNumThread()).wait();
        return temp.mean();
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class ForceModel>
    T RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::calcPotentialClassical(const ForceModel& model) const {
        T result = 0;
        for (size_t i = 0; i < getNumReplica(); ++i)
            result += model.potentialV(phaseToCell(i));
        return result;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    T RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::calcClassicalElastic() const {
        const size_t dof = getDOF();
        auto pos = getPhaseMatrix().bottomRows(dof);
        const T omegaW = ringPolymer.calcOmegaW(temperatureT);
        T result = 0;
        for (size_t i = 0; i < dof; ++i) {
            const T mass = cell.getMass(i / Dim);
            for (size_t j = 0; j < getNumReplica(); ++j)
                result += mass * square(omegaW * (pos(i, j) - pos(i, (j + 1) % getNumReplica()))) * 0.5;
        }
        return result;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class ForceModel>
    T RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::calcClassicalInternalEnergy(const ForceModel& model) const {
        return calcKineticClassical() + calcPotentialClassical<ForceModel>(model) + calcClassicalElastic();
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel>
    T RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::calcTemperature() const {
        return ringPolymer.template calcTemperature<KineticModel>();
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel, class ForceModel, class Executor>
    T RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::calcPressThermo(ForceModel& model) const {
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        constexpr bool IsPeriodBoundary = Traits<ForceModel>::IsPeriodBoundary;
        T result = calcKinetic<KineticModel>() / (getVolume() * (Dim * 0.5));
        if constexpr (!isFreeModel) {
            const size_t numReplica = getNumReplica();
            T temp = 0;
            for (size_t i = 0; i < numReplica; ++i) {
                auto cell = phaseToCell(i);
                if constexpr (IsPeriodBoundary)
                    cell.normalize();
                temp += model.virial(std::move(cell)).trace();
            }
            result += temp / T(numReplica * Dim);
        }
        return result;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class ForceModel, class Executor>
    RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::LatticeMatrix
    RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::makeStressPrim(ForceModel& model) const {
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        if constexpr (NumReplica == 1)
            return makeStressClassical<ForceModel, Executor>(model);

        Array<LatticeMatrix> buffer(getNumReplica());
        const T squaredOmegaW = square(ringPolymer.calcOmegaW(temperatureT));
        auto kernel = [this, &model, &buffer, squaredOmegaW](unsigned int replica) {
            using VectorType = DenseVector<T, Dim>;
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
                const T mass = getMassVec()[i];
                const T repMass = reciprocal(mass);
                const auto atomMomentum = momentum.segment(from, to);
                kineticStress += (repMass * atomMomentum) * atomMomentum.transpose();

                const auto atomPos = pos.segment(from, to);
                const auto atomPos1 = pos1.segment(from, to);
                const VectorType deltaPos = atomPos - atomPos1;
                const T factorK = mass * squaredOmegaW;
                kineticStress -= (factorK * deltaPos) * deltaPos.transpose();
            }
            buffer[replica] = kineticStress * reciprocal(getVolume());

            std::ignore = model; // Silent warnings
            if constexpr (!isFreeModel) {
                constexpr bool IsPeriodBoundary = Traits<ForceModel>::IsPeriodBoundary;
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
    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class ForceModel, class Executor>
    RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::LatticeMatrix
    RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::makeStressVirial(ForceModel& model) const {
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        static_assert(!isFreeModel, "[Error]: This function does not apply to ideal gas model");
        if constexpr (NumReplica == 1)
            return makeStressClassical<ForceModel, Executor>(model);

        Array<LatticeMatrix> buffer(getNumReplica());
        const auto centroidPos = ringPolymer.makeCentroidPos();
        auto kernel = [this, &model, &buffer, &centroidPos](unsigned int replica) {
            using VectorType = DenseVector<T, Dim>;
            const size_t dof = getDOF();
            const size_t numReplica = getNumReplica();
            const auto col = getPhaseMatrix().col(replica);
            const auto momentum = col.head(dof);
            const auto pos = col.tail(dof);
            const auto centroid = centroidPos.flatten();
            const auto force = forceBuffer.col(replica);

            LatticeMatrix kineticStress(Dim, Dim, 0);
            for (size_t i = 0; i < getNumParticle(); ++i) {
                const size_t from = i * Dim;
                const size_t to = from + Dim;
                const T mass = getMassVec()[i];
                const T factor = reciprocal(mass * T(numReplica));
                const auto atomMomentum = momentum.segment(from, to);
                kineticStress += (factor * atomMomentum) * atomMomentum.transpose();

                const auto atomPos = pos.segment(from, to);
                const auto atomCentroid = centroid.segment(from, to);
                const auto atomForce = force.segment(from, to);
                const VectorType deltaPos = atomPos - atomCentroid;
                kineticStress -= deltaPos * atomForce.transpose();
            }
            buffer[replica] = kineticStress * reciprocal(getVolume());

            std::ignore = model; // Silent warnings
            if constexpr (!isFreeModel) {
                constexpr bool IsPeriodBoundary = Traits<ForceModel>::IsPeriodBoundary;
                auto cell = phaseToCell(replica);
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

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class ForceModel, class Executor>
    RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::LatticeMatrix
    RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::makeStressClassical(ForceModel& model) const {
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        constexpr bool IsPeriodBoundary = Traits<ForceModel>::IsPeriodBoundary;
        LatticeMatrix result(Dim, Dim, 0);
        if constexpr (NumReplica == 1) {
            const T repVolume = reciprocal(getVolume());
            const auto col = getPhaseMatrix().col(0);
            const auto momentum = col.head(getDOF());
            for (size_t i = 0; i < getNumParticle(); ++i) {
                const size_t from = i * Dim;
                const size_t to = from + Dim;
                const T repMass = reciprocal(getMassVec()[i]);
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
            Array<LatticeMatrix> buffer(getNumReplica());
            auto kernel = [this, &model, &buffer](unsigned int replica) {
                const size_t dof = getDOF();
                const size_t numReplica = getNumReplica();
                const auto col = getPhaseMatrix().col(replica);
                const auto momentum = col.head(dof);

                LatticeMatrix kineticStress(Dim, Dim, 0);
                for (size_t i = 0; i < getNumParticle(); ++i) {
                    const size_t from = i * Dim;
                    const size_t to = from + Dim;
                    const T mass = getMassVec()[i];
                    const T factor = reciprocal(mass * T(numReplica));
                    const auto atomMomentum = momentum.segment(from, to);
                    kineticStress += (factor * atomMomentum) * atomMomentum.transpose();
                }
                buffer[replica] = kineticStress * reciprocal(getVolume());

                std::ignore = model; // Silent warnings
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

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel, class ForceModel, class Executor>
    VectorND<T> RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::testNVE(
            T duration, KineticModel& kineticModel, ForceModel& forceModel) const {
        This rpmd = *this;
        const uint64_t step = Base::durationToStep(duration, timeStep);
        VectorND<T> pot(step);
        for (uint64_t i = 0; i < step; ++i) {
            rpmd.nve_step<KineticModel, ForceModel, Executor>(kineticModel, forceModel);
            pot[i] = rpmd.calcKinetic<KineticModel>() + rpmd.calcPotential<ForceModel, Executor>(forceModel);
        }
        pot -= T(pot[0]);
        return pot;
    }
#ifdef PHYSICA_HDF5
    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::read(const H5Loc& loc, const char* name) {
        const auto group = loc.openGroup(name);
        LatticeMatrix temp{};
        temp.read(group, "lattice");
        setLattice(std::move(temp));
        getPhaseMatrix().read(group, "phase");
        forceBuffer.read(group, "force");
    }
    
    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::write(H5Loc& loc, const char* name) const {
        auto group = loc.openGroup(name);
        getLattice().write(group, "lattice");
        getPhaseMatrix().write(group, "phase");
        forceBuffer.write(group, "force");
    }
#endif
    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::swap(RPMD& __restrict obj) noexcept {
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

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    inline void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::setTemperature(T temperature) {
        assert(!temperature.isNegative() && "[Error]: Negative temperature is not physical");
        temperatureT = temperature;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::setTimeStep(T timeStep_) {
        timeStep = timeStep_;
        checkParam();
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::toContractBeadRepr(size_t posID) {
        assert(posID < getDOF());
        auto row = ringPolymer.getBuffer().row(1);
        auto head = row.head(fftContract.getKSpaceSize());
        fftContract.invTransform(head);
        auto pos = posContract.row(posID);
        pos = fftContract.getRSpace() * (T(getNumContract()) / T(getNumReplica()));
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::forceToNormRepr(size_t posID) {
        assert(posID < getDOF());
        fftContract.transform(forceContract.row(posID));
        auto row = ringPolymer.getBuffer().row(0);
        row = T(0);
        auto head = row.head(fftContract.getKSpaceSize());
        head = fftContract.getKSpace();
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::forceToBeadRepr(size_t posID) {
        assert(posID < getDOF());
        auto& buffer = ringPolymer.getBuffer();
        auto& fft = ringPolymer.getCanonicalFFT();
        fft.invTransform(buffer.row(0));
        auto f = forceBuffer.row(posID);
        f += fft.getRSpace() * (T(getNumReplica()) / T(getNumContract()));
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::contract() {
        const size_t dof = getDOF();
        for (size_t i = 0; i < dof; ++i) {
            ringPolymer.toNormalRepr(i);
            toContractBeadRepr(i);
        }
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::decontract() {
        const size_t dof = getDOF();
        for (size_t i = 0; i < dof; ++i) {
            forceToNormRepr(i);
            forceToBeadRepr(i);
        }
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::forceStep(T deltaT) {
        const size_t dof = getDOF();
        for (size_t replica = 0; replica < getNumReplica(); ++replica) {
            auto col = getPhaseMatrix().col(replica);
            auto momentum = col.head(dof);
            momentum += forceBuffer.col(replica) * deltaT;
        }
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    bool RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::checkCentroid() const {
        constexpr bool isGood = true;
        PositionMatrix centroid = ringPolymer.makeCentroidPos();
        cell.toDirect(centroid);
        for (auto& elem : centroid)
            if (!(T(0) <= elem && elem <= T(0)))
                return !isGood;
        return isGood;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::checkParam() const {
        if constexpr (NumReplica != 1) {
            const T cycle = Tv(2 * M_PI) / ringPolymer.calcOmegaW(temperatureT);
            bool isSmallEnough = timeStep < cycle / Tv(4);
            if (!isSmallEnough)
                throw std::invalid_argument("[Error]: Time step is too large");
        }
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    uint64_t RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::durationToStep(Tv duration, Tv timeStep) {
        return double(duration / timeStep) + 0.5;
    }
}
