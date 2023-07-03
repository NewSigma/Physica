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

#include "Physica/Core/Math/Transform/FFT.h"
#include "Physica/Core/Math/Statistics/NumCharacter.h"
#include "Physica/Core/Physics/PhyConst.h"
#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"
#include "MDCell.h"
#include "MDImpl/RingPolymer.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] Habershon S, Manolopoulos D E, Markland T E, et al. Ring-Polymer Molecular Dynamics: Quantum Effects in Chemical Dynamics from Classical Trajectories in an Extended Phase Space[J]. Annual Review of Physical Chemistry, 2013, 64(1):387-413.
     * [2] Jos Thijssen. Computational Physics[M].London: Cambridge university press, 2013:197-211
     * [3] Liu J, Li D, Liu X. A simple and accurate algorithm for path integral molecular dynamics with the Langevin thermostat[J]. J. Chem. Phys, 2016, 145(2):1291-1301.
     * [4] T. E. Markland, D. E. Manolopoulos. An efficient ring polymer contraction scheme for imaginary time path integral simulations[J]. J. Chem. Phys. 129, 024105 (2008)
     * [5] Ian R. Craig and David E. Manolopoulos, J. Chem. Phys. 121, 3368 (2004)
     * 
     * TODO: replace several ScalarType to PosScalarType
     */
    template<class ScalarType, class PosScalarType, unsigned int Dim = 3, size_t NumReplica = Dynamic>
    class RPMD final {
    public:
        using RingPolymerType = RingPolymer<ScalarType, PosScalarType, Dim, NumReplica>;
        using PhaseMatrix = typename RingPolymerType::PhaseMatrix;
        using ForceMatrix = DenseMatrix<ScalarType, MatrixOption::Column | MatrixOption::Vector>;
        using MDCellType = typename RingPolymerType::MDCellType;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using PositionMatrix = typename MDCellType::PositionMatrix;
    private:
        MDCellType cell;
        RingPolymerType ringPolymer;
        ForceMatrix forceBuffer;

        FFT<PosScalarType, 1> fftContract;
        PhaseMatrix posContract;
        ForceMatrix forceContract;
        /* Constant */
        ScalarType temperatureT;
        ScalarType timeStep;
    public:
        RPMD(MDCellType cell_,
             size_t numReplica,
             size_t numContract,
             ScalarType temperatureT_,
             ScalarType timeStep_);
        RPMD(const RPMD&) = default;
        RPMD(RPMD&&) noexcept = default;
        ~RPMD() = default;
        /* Operators */
        RPMD& operator=(RPMD obj) noexcept;
        /* Operations */
        template<class ForceModel, class Executor> void updateForce(const ForceModel& model);
        template<class KineticModel,
                 class ForceModel,
                 class Executor>
        void nve_step(KineticModel& kineticModel, const ForceModel& forceModel);
        template<class KineticModel,
                 class ForceModel,
                 class Executor>
        void nve_step_for(ScalarType duration, KineticModel& kineticModel, const ForceModel& forceModel);

        template<class Thermostat,
                 class RandomGenerator,
                 class KineticModel,
                 class ForceModel,
                 class Executor>
        void nvt_step(const Thermostat& thermostat, RandomGenerator& gen, KineticModel& kineticModel, const ForceModel& forceModel);
        template<class Thermostat,
                 class RandomGenerator,
                 class KineticModel,
                 class ForceModel,
                 class Executor>
        void nvt_step_for(ScalarType duration, const Thermostat& thermostat, RandomGenerator& gen, KineticModel& kineticModel, const ForceModel& forceModel);

        template<class Thermostat,
                 class RandomGenerator,
                 class Barostat,
                 class KineticModel,
                 class ForceModel,
                 class Executor>
        void npt_step(const Thermostat& thermostat, RandomGenerator& gen, Barostat& barostat, KineticModel& kineticModel, const ForceModel& forceModel);
        template<class Thermostat,
                 class RandomGenerator,
                 class Barostat,
                 class KineticModel,
                 class ForceModel,
                 class Executor>
        void npt_step_for(ScalarType duration, const Thermostat& thermostat, RandomGenerator& gen, Barostat& barostat, KineticModel& kineticModel, const ForceModel& forceModel);

        template<class RandomGenerator> void initMomentum(RandomGenerator& gen);
        void scaleVelocity();
        void normalizeCentroid();
        [[nodiscard]] MDCellType phaseToCell(size_t replica) const;
        [[nodiscard]] MDCellType contractToCell(size_t contract) const;
        [[nodiscard]] MDCellType makeAverageCell() const;
        void checkParam() const;
        void swap(RPMD& obj) noexcept;
        /* Getters */
        [[nodiscard]] constexpr unsigned int getDim() const noexcept { return Dim; }
        [[nodiscard]] const LatticeMatrix& getLattice() const noexcept { return cell.getLattice(); }
        [[nodiscard]] const typename MDCellType::InvLatticeMatrix& getInvLattice() const noexcept { return cell.getInvLattice(); }
        [[nodiscard]] const typename MDCellType::MassVector& getMassVec() const noexcept { return cell.getMassVec(); }
        [[nodiscard]] size_t getNumParticle() const noexcept { return cell.getNumParticle(); }
        [[nodiscard]] PosScalarType getVolume() const noexcept { return cell.getVolume(); }
        [[nodiscard]] const RingPolymerType& getRingPolymer() const noexcept { return ringPolymer; }
        [[nodiscard]] RingPolymerType& getRingPolymer() noexcept { return ringPolymer; }
        [[nodiscard]] const PhaseMatrix& getPhaseMatrix() const noexcept { return ringPolymer.asMatrix(); }
        [[nodiscard]] PhaseMatrix& getPhaseMatrix() noexcept { return ringPolymer.asMatrix(); }
        [[nodiscard]] size_t getDOF() const noexcept { return ringPolymer.getDOF(); }
        [[nodiscard]] size_t getNumReplica() const noexcept { return ringPolymer.getNumReplica(); }
        [[nodiscard]] size_t getKSpaceSize() const noexcept { return ringPolymer.getKSpaceSize(); }
        [[nodiscard]] size_t getNumContract() const noexcept { return fftContract.getRSpaceSize(); }
        [[nodiscard]] bool isContractEnabled() const noexcept { return getNumReplica() != getNumContract(); }
        [[nodiscard]] ScalarType getTemperature() const noexcept { return temperatureT; }

        [[nodiscard]] const ForceMatrix& getForce() const noexcept { return forceBuffer; }

        template<class ForceModel>
        [[nodiscard]] ScalarType getClassicalPotentialEnergy(const ForceModel& model) const;
        [[nodiscard]] ScalarType getClassicalElastic() const;
        template<class ForceModel>
        [[nodiscard]] ScalarType getClassicalInternalEnergy(const ForceModel& model) const;
        [[nodiscard]] ScalarType calcKinetic() const;
        template<class ForceModel, class Executor>
        [[nodiscard]] ScalarType calcPotential(const ForceModel& model) const;

        template<class ForceModel>
        [[nodiscard]] LatticeMatrix makeStress(const ForceModel& model) const;
        /* Setters */
        void setTemperature(ScalarType temperature);
    private:
        void toContractBeadRepr(size_t posID);
        void forceToNormRepr(size_t posID);
        void forceToBeadRepr(size_t posID);
        void contract();
        void decontract();
        void forceStep(ScalarType deltaT);
        bool checkCentroid() const;
    };
}

#include "MDImpl/RPMDImpl.h"
