/*
 * Copyright 2022-2026 Weibo He.
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
            , fftContract(numContract, PlanFlag::Estimate) {
        assert(0 < numContract && numContract <= numReplica);
        assert(NumReplica == Dynamic || NumReplica == numReplica);
        ringPolymer = RingPolymerType(cell, numReplica);

        const size_t dof = getDOF();
        forceBuffer.resize(dof, numReplica, 0);
        if constexpr (!IsClassical) {
            if (isContractEnabled()) {
                posContract.resize(dof, numContract);
                forceContract.resize(dof, numContract);
            }
        }
        setTemperature(temperatureT_);
        setTimeStep(timeStep_);
    }
    /**
     * Contract method to improve performance introduced in [1]
     * 
     * Reference:
     * [1] J. Chem. Phys. 129, 024105 (2008); https://doi.org/10.1063/1.2953308
     */
    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<ExecutePolicy P>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::updateForce(auto& model) {
        using ForceModel = std::remove_cvref_t<decltype(model)>;
        constexpr bool IsPeriodBoundary = model.isPeriodBoundary();
        // TODO: Stateful ForceModel is not parallelizable. Add a trait to determine whether it is stateful.
        // Currently we assume GPU models are stateful
        // constexpr ExecutePolicy HostPolicy = Traits<ForceModel>::IsStateful ? Sequential : Thread;
        constexpr ExecutePolicy HostPolicy = P == GPU ? Sequential : Thread;
        static_assert((P != GPU) || std::allocator_traits<ForceMatrixAllocator>::isPageLocked, "[Error]: Allocator is not page locked, performance will decrease");
        if (!isContractEnabled()) {
            auto kernel = [this, &model](unsigned int replica) {
                MDCellType cell = phaseToCell(replica);
                if constexpr (IsPeriodBoundary)
                    cell.normalize();
                auto saveTo = forceBuffer.col(replica);
                model.template forceAsync<P>(std::move(cell), saveTo);
            };

            if constexpr (NumReplica == 1)
                kernel(0);
            else
                parallel_for<HostPolicy>(kernel, getNumReplica()).wait();

            if constexpr (P == GPU)
                Task<P>::wait();
            return;
        }

        if constexpr (!IsClassical) {
            if constexpr (Traits<ForceModel>::IsContractable) {
                auto task_uncontract = parallel_for<HostPolicy>([&](size_t replica) {
                    MDCellType cell = phaseToCell(replica);
                    if constexpr (IsPeriodBoundary)
                        cell.normalize();
                    auto saveTo = forceBuffer.col(replica);
                    saveTo = model.template force_uncontract<P>(std::move(cell));
                }, getNumReplica(), 0);

                contract();
                auto task_contract = parallel_for<HostPolicy>([&](size_t contract) {
                    MDCellType cell = contractToCell(contract);
                    if constexpr (IsPeriodBoundary)
                        cell.normalize();
                    auto saveTo = forceContract.col(contract);
                    saveTo = model.template force_contract<P>(std::move(cell));
                }, getNumContract(), 0);

                task_uncontract.wait();
                task_contract.wait();
                decontract();
            }
            else
                noImpl("[Error]: Force contract is not implemented");
        }
        else
            unreachable();
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<ExecutePolicy P>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::nve_step(auto& kineticModel, auto& forceModel) {
        nve_step_impl<P>(kineticModel, forceModel, timeStep);
    }
    /**
     * Some applications may need to step forward and back, e.g. \class HamiltonMC
     */
    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<ExecutePolicy P>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::nve_step_back(auto& kineticModel, auto& forceModel) {
        nve_step_impl<P>(kineticModel, forceModel, -timeStep);
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<ExecutePolicy P>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::nve_step_for(
            T duration, auto& kineticModel, auto& forceModel) {
        const uint64_t step = durationToStep(duration, timeStep);
        for (uint64_t _ = 0; _ < step; ++_)
            nve_step<P>(kineticModel, forceModel);
    }
    /**
     * BAOAB integrator as introduced in [1]
     * 
     * Reference:
     * [1] J. Chem. Phys, 145, 024103 (2016); https://doi.org/10.1063/1.4954990
     */
    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<RNG R, ExecutePolicy P>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::nvt_step(const auto& thermostat, auto& kineticModel, auto& forceModel) {
        using Thermostat = std::remove_cvref_t<decltype(thermostat)>;
        using ForceModel = std::remove_cvref_t<decltype(forceModel)>;
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        constexpr bool IsPeriodBoundary1 = kineticModel.isPeriodBoundary();
        constexpr bool IsPeriodBoundary2 = forceModel.isPeriodBoundary();
        static_assert(isFreeModel || (IsPeriodBoundary1 == IsPeriodBoundary2), "[Error]: Inconsistent boundary condition");

        constexpr bool isSeedFixed = Traits<R>::IsSeedFixed;
        constexpr auto NoRandPolicy = isSeedFixed ? Sequential : P;

        constexpr bool IsCentroidCoupled = Traits<Thermostat>::IsCentroidCoupled;
        if constexpr (IsCentroidCoupled && IsPeriodBoundary1)
            assert(thermostat.isRemoveDriftEnabled() && "[Error]: Because the KineticModel has period boundary, thermostat should remove its effect on centroid");

        if constexpr (isFreeModel) {
            kineticModel.nve_step(ringPolymer, timeStep * 0.5);
            thermostat.template step<R, NoRandPolicy>(ringPolymer, timeStep);
            kineticModel.nve_step(ringPolymer, timeStep * 0.5);
        }
        else {
            forceStep(timeStep * 0.5);
            kineticModel.nve_step(ringPolymer, timeStep * 0.5);
            thermostat.template step<R, NoRandPolicy>(ringPolymer, timeStep);
            kineticModel.nve_step(ringPolymer, timeStep * 0.5);
            updateForce<P>(forceModel);
            forceStep(timeStep * 0.5);
        }
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<RNG R, ExecutePolicy P>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::nvt_step_for(T duration, const auto& thermostat, auto& kineticModel, auto& forceModel) {
        const uint64_t step = durationToStep(duration, timeStep);
        for (uint64_t _ = 0; _ < step; ++_)
            nvt_step<R, P>(thermostat, kineticModel, forceModel);
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<RNG R, ExecutePolicy P>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::npt_step(const auto& thermostat, auto& barostat, auto& kineticModel, auto& forceModel) {
        using Thermostat = std::remove_cvref_t<decltype(thermostat)>;
        using Barostat = std::remove_cvref_t<decltype(barostat)>;
        using ForceModel = std::remove_cvref_t<decltype(forceModel)>;
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        constexpr bool IsPeriodBoundary1 = kineticModel.isPeriodBoundary();
        constexpr bool IsPeriodBoundary2 = forceModel.isPeriodBoundary();
        static_assert(isFreeModel || (IsPeriodBoundary1 == IsPeriodBoundary2), "[Error]: Inconsistent boundary condition");

        constexpr unsigned int BarostatOrder = Traits<Barostat>::Order;
        static_assert(BarostatOrder == 2 || BarostatOrder == 1, "[Error]: Invalid barostat");

        constexpr bool isSeedFixed = Traits<R>::IsSeedFixed;
        constexpr auto NoRandPolicy = isSeedFixed ? Sequential : P;

        constexpr bool IsCentroidCoupled = Traits<Thermostat>::IsCentroidCoupled;
        if constexpr (IsCentroidCoupled && IsPeriodBoundary1)
            assert(thermostat.isRemoveDriftEnabled() && "[Error]: Because the KineticModel has period boundary, thermostat should remove its effect on centroid");

        if constexpr (BarostatOrder == 2) {
            barostat.forceStep(*this, forceModel, timeStep * 0.5);
            kineticModel.npt_step(*this, barostat, timeStep * 0.5);
            thermostat.template step<R, NoRandPolicy>(ringPolymer, timeStep);
            kineticModel.npt_step(*this, barostat, timeStep * 0.5);
            updateForce<P>(forceModel);
            barostat.forceStep(*this, forceModel, timeStep * 0.5);
        }
        else {
            const LatticeMatrix stress = makeStressPrim<P>(forceModel);
            forceStep(timeStep * 0.5);
            kineticModel.nve_step(ringPolymer, timeStep * 0.5);
            barostat.template npt_step<ForceModel>(*this, stress, timeStep * 0.5);
            thermostat.template step<R, NoRandPolicy>(ringPolymer, timeStep);
            barostat.template npt_step<ForceModel>(*this, stress, timeStep * 0.5);
            kineticModel.nve_step(ringPolymer, timeStep * 0.5);
            if constexpr (Traits<ForceModel>::IsLatticeDependent)
                forceModel.setLattice(getLattice());
            updateForce<P>(forceModel);
            forceStep(timeStep * 0.5);
        }
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<RNG R, ExecutePolicy P>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::npt_step_for(T duration, const auto& thermostat, auto& barostat, auto& kineticModel, auto& forceModel) {
        const uint64_t step = durationToStep(duration, timeStep);
        for (uint64_t _ = 0; _ < step; ++_)
            npt_step<R, P>(thermostat, barostat, kineticModel, forceModel);
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
    template<ExecutePolicy P>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::fire_vstep(FireModel<T, Dim>& fire, auto& kineticModel, auto& forceModel) {
        using KineticModel = std::remove_cvref_t<decltype(kineticModel)>;
        using ForceModel = std::remove_cvref_t<decltype(forceModel)>;
        checkRelaxParam<KineticModel, ForceModel>();

        kineticModel.nve_step(ringPolymer, timeStep);
        updateForce<P>(forceModel);
        fire.paramStep(*this);
        forceStep(fire.getTimeStep());
        fire.mixingStep(*this);
    }
    /**
     * fire_pstep is fire_p(ress)step
     */
    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<BaroType Type, ExecutePolicy P>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::fire_pstep(CFireModel<T, Dim, Type>& cfire, auto& kineticModel, auto& forceModel) {
        using KineticModel = std::remove_cvref_t<decltype(kineticModel)>;
        using ForceModel = std::remove_cvref_t<decltype(forceModel)>;
        checkRelaxParam<KineticModel, ForceModel>();

        kineticModel.nve_step(ringPolymer, timeStep);
        cfire.nve_step(*this);
        if constexpr (Traits<ForceModel>::IsLatticeDependent)
            forceModel.setLattice(getLattice());
        updateForce<P>(forceModel);
        const LatticeMatrix stress = forceModel.virial(phaseToCell(0));
        cfire.paramStep(*this);
        forceStep(cfire.getTimeStep());
        cfire.forceStep(*this, std::move(stress));
        cfire.mixingStep(*this);
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel, RNG R>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::initMomentum() {
        return ringPolymer.template initMomentum<KineticModel, R>(temperatureT);
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::scaleVelocity() {
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
    auto RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::phaseToCell(size_t replica) const -> MDCellType {
        return MDCellType(cell.getLattice(), ringPolymer.makeBeadPos(replica), cell.getMassVec());
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    auto RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::contractToCell(size_t contract) const -> MDCellType {
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
    auto RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::makeAverageCell() const -> MDCellType {
        return MDCellType(getLattice(), ringPolymer.makeCentroidPos(), cell.getMassVec());
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel>
    T RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::calcKinetic() const {
        if constexpr (IsClassical)
            return calcKineticClassical();
        else {
            /**
             * Kinetic energy using virial estimator referenced from [1]. Generally, we should use this one.
             * 
             * Reference:
             * [1] J. Chem. Phys. 76, 5150-5155 (1982); https://doi.org/10.1063/1.442815
             */
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
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel>
    T RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::calcKinetic(size_t dofIndex) const {
        static_assert(!IsClassical, "[Error]: Use classical kinetic instead");
        const T repBeta = ringPolymer.calcRepBeta(calcTemperature<KineticModel>());
        const auto pos = getPhaseMatrix().row(getDOF() + dofIndex);
        const T centroidPos = pos.mean();

        T kinetic = repBeta - (pos - centroidPos) * forceBuffer.row(dofIndex);
        kinetic /= T(getNumReplica() * 2);
        return kinetic;
    }
    /**
     * Kinetic energy using primitive estimator referenced from [1]
     *
     * Note: Prefer this estimator either
     * 1. NumReplica is small -- prim estimator has lower variance
     * 2. force model is \class EmptyForceModel -- virial estimator does not work
     * 
     * Reference:
     * [1] J. Chem. Phys. 76, 5150 (1982); https://doi.org/10.1063/1.442815
     */
    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel>
    T RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::calcKineticPrim() const {
        static_assert(!IsClassical, "[Error]: Use classical kinetic instead");
        T kinetic = 0;
        for (size_t i = 0; i < getDOF(); ++i)
            kinetic += calcKineticPrim<KineticModel>(i);
        return kinetic;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel>
    T RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::calcKineticPrim(size_t dofIndex) const {
        static_assert(!IsClassical, "[Error]: Use classical kinetic instead");
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
    T RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::calcKineticClassical() const {
        return ringPolymer.calcKineticClassical();
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<ExecutePolicy P>
    T RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::calcPotential(auto& forceModel) const {
        VectorND<T> temp(getNumReplica());
        auto kernel = [this, &forceModel, &temp](unsigned int replica) {
            temp[replica] = forceModel.potentialV(phaseToCell(replica));
        };
        parallel_for<P>(kernel, getNumReplica(), 0).wait();
        return temp.mean();
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<ExecutePolicy P>
    T RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::calcPotentialClassical(auto& forceModel) const {
        T result = 0;
        for (size_t i = 0; i < getNumReplica(); ++i)
            result += forceModel.template potentialV<P>(phaseToCell(i));
        return result;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    T RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::calcClassicalElastic() const {
        if constexpr (IsClassical)
            return 0;
        else {
            T result = 0;
            const size_t dof = getDOF();
            auto pos = getPhaseMatrix().bottomRows(dof);
            const T omegaW = ringPolymer.calcOmegaW(temperatureT);
            for (size_t i = 0; i < dof; ++i) {
                const T mass = cell.getMass(i / Dim);
                for (size_t j = 0; j < getNumReplica(); ++j)
                    result += mass * square(omegaW * (pos[i, j] - pos[i, (j + 1) % getNumReplica()])) * 0.5;
            }
            return result;
        }
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<ExecutePolicy P>
    T RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::calcClassicalInternalEnergy(auto& forceModel) const {
        return calcKineticClassical() + calcPotentialClassical<P>(forceModel) + calcClassicalElastic();
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel>
    T RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::calcTemperature() const {
        return ringPolymer.template calcTemperature<KineticModel>();
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<class KineticModel, ExecutePolicy P>
    T RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::calcPressThermo(auto& forceModel) const {
        using ForceModel = std::remove_cvref_t<decltype(forceModel)>;
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        constexpr bool IsPeriodBoundary = forceModel.isPeriodBoundary();
        T result = calcKinetic<KineticModel>() / (getVolume() * (Dim * 0.5));
        if constexpr (!isFreeModel) {
            const size_t numReplica = getNumReplica();
            T temp = 0;
            for (size_t i = 0; i < numReplica; ++i) {
                auto cell = phaseToCell(i);
                if constexpr (IsPeriodBoundary)
                    cell.normalize();
                temp += forceModel.virial(std::move(cell)).trace();
            }
            result += temp / T(numReplica * Dim);
        }
        return result;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<ExecutePolicy P>
    auto RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::makeStressPrim(auto& forceModel) const -> LatticeMatrix {
        using ForceModel = std::remove_cvref_t<decltype(forceModel)>;
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        if constexpr (NumReplica == 1)
            return makeStressClassical<P>(forceModel);

        Array<LatticeMatrix> buffer(getNumReplica());
        const T squaredOmegaW = square(ringPolymer.calcOmegaW(temperatureT));
        auto kernel = [this, &forceModel, &buffer, squaredOmegaW](unsigned int replica) {
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

            std::ignore = forceModel; // Silent warnings
            if constexpr (!isFreeModel) {
                constexpr bool IsPeriodBoundary = forceModel.isPeriodBoundary();
                MDCellType cell = phaseToCell(replica);
                if constexpr (IsPeriodBoundary)
                    cell.normalize();
                buffer[replica] += forceModel.virial(cell);
            }
        };
        parallel_for<P>(kernel, getNumReplica(), 0).wait();
        LatticeMatrix result(Dim, Dim, 0);
        for (size_t i = 0; i < buffer.getLength(); ++i)
            result.toNextMean(i, buffer[i]);
        return result;
    }
    /**
     * Reference:
     * [1] Comp. Phys. Comm. 185, 1019 (2013); https://doi.org/10.1016/j.cpc.2013.10.027
     */
    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<ExecutePolicy P>
    auto RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::makeStressVirial(auto& forceModel) const -> LatticeMatrix {
        using ForceModel = std::remove_cvref_t<decltype(forceModel)>;
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        static_assert(!isFreeModel, "[Error]: This function does not apply to ideal gas model");
        if constexpr (NumReplica == 1)
            return makeStressClassical<P>(forceModel);

        Array<LatticeMatrix> buffer(getNumReplica());
        const auto centroidPos = ringPolymer.makeCentroidPos();
        auto kernel = [this, &forceModel, &buffer, &centroidPos](unsigned int replica) {
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

            std::ignore = forceModel; // Silent warnings
            if constexpr (!isFreeModel) {
                constexpr bool IsPeriodBoundary = forceModel.isPeriodBoundary();
                auto cell = phaseToCell(replica);
                if constexpr (IsPeriodBoundary)
                    cell.normalize();
                buffer[replica] += forceModel.virial(cell);
            }
        };
        parallel_for<P>(kernel, getNumReplica(), 0).wait();
        LatticeMatrix result(Dim, Dim, 0);
        for (size_t i = 0; i < buffer.getLength(); ++i)
            result.toNextMean(i, buffer[i]);
        return result;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<ExecutePolicy P>
    auto RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::makeStressClassical(auto& forceModel) const -> LatticeMatrix {
        using ForceModel = std::remove_cvref_t<decltype(forceModel)>;
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        constexpr bool IsPeriodBoundary = forceModel.isPeriodBoundary();
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
                result += forceModel.virial(cell);
            }
        }
        else {
            Array<LatticeMatrix> buffer(getNumReplica());
            auto kernel = [this, &forceModel, &buffer](unsigned int replica) {
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

                std::ignore = forceModel; // Silent warnings
                if constexpr (!isFreeModel) {
                    MDCellType cell = phaseToCell(replica);
                    if constexpr (IsPeriodBoundary)
                        cell.normalize();
                    buffer[replica] += forceModel.virial(cell);
                }
            };
            parallel_for<P>(kernel, getNumReplica(), 0).wait();
            for (size_t i = 0; i < buffer.getLength(); ++i)
                result.toNextMean(i, buffer[i]);
        }
        return result;
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
    auto&& RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::getRingPolymer(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).ringPolymer;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    auto&& RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::getPhaseMatrix(this auto&& self) noexcept {
        return self.getRingPolymer().asMatrix();
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    bool RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::isContractEnabled() const noexcept {
        if constexpr (IsClassical)
            return false;
        else
            return getNumReplica() != getNumContract();
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::setTemperature(T temperature) noexcept {
        assert(!temperature.isNegative() && "[Error]: Negative temperature is not physical");
        temperatureT = temperature;
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::setTimeStep(T timeStep_) {
        assert(timeStep_.isPositive());
        timeStep = timeStep_;
        if constexpr (!IsClassical) {
            const T cycle = Tv(2 * M_PI) / ringPolymer.calcOmegaW(temperatureT);
            bool isSmallEnough = timeStep < cycle / Tv(4);
            if (!isSmallEnough)
                throw std::invalid_argument("[Error]: Time step is too large");
        }
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    uint64_t RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::durationToStep(Tv duration, Tv timeStep) {
        return std::lround((duration / timeStep).toMachine());
    }

    template<Scalar T, unsigned int Dim, size_t NumReplica, class ForceMatrixAllocator>
    template<ExecutePolicy P>
    void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::nve_step_impl(auto& kineticModel, auto& forceModel, Tv deltaT) {
        using ForceModel = std::remove_cvref_t<decltype(forceModel)>;
        constexpr bool isFreeModel = Internal::is_empty_force_model<ForceModel>::value;
        constexpr bool IsPeriodBoundary1 = kineticModel.isPeriodBoundary();
        constexpr bool IsPeriodBoundary2 = forceModel.isPeriodBoundary();
        static_assert(isFreeModel || (IsPeriodBoundary1 == IsPeriodBoundary2), "[Error]: Inconsistent boundary condition");
        if (isFreeModel)
            kineticModel.nve_step(ringPolymer, deltaT);
        else {
            forceStep(deltaT * Tv(0.5));
            kineticModel.nve_step(ringPolymer, deltaT);
            updateForce<P>(forceModel);
            forceStep(deltaT * Tv(0.5));
        }
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
    template<class KineticModel, class ForceModel>
    constexpr void RPMD<T, Dim, NumReplica, ForceMatrixAllocator>::checkRelaxParam() {
        constexpr bool IsPeriodBoundary1 = KineticModel::isPeriodBoundary();
        constexpr bool IsPeriodBoundary2 = ForceModel::isPeriodBoundary();
        static_assert(IsPeriodBoundary1 == IsPeriodBoundary2, "[Error]: Inconsistent boundary condition");
        static_assert(NumReplica == 1, "[Error]: Relaxing using PIMD makes no sence, NumReplica = 1 shall be enough");
        static_assert(!Internal::is_empty_force_model<ForceModel>::value, "[Error]: Relax a empty model does nothing");
    }
}
