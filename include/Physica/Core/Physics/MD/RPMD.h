/*
 * Copyright 2022-2024 Weibo He.
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

#include <Physica/Core/Math/Transform/FFT.h>
#include <Physica/Core/Math/Statistics/NumCharacter.h>
#include <Physica/Core/Physics/PhyConst.h>
#include <Physica/Core/Parallel/Executor/SequentialExecutor.h>
#include "MDCell.h"
#include "MDImpl/RingPolymer.h"
#include "Barostat/BaroType.h"

namespace Physica::Core {
    template<class ScalarType, unsigned int Dim> class FireModel;
    template<class ScalarType, unsigned int Dim, BaroType Type> class CFireModel;

    template<class ScalarType>
    class RPMDBase {        
    public:
        using PlainScalar = typename ScalarType::PlainScalar;

        [[nodiscard]] static uint64_t durationToStep(PlainScalar duration, PlainScalar timeStep) {
            return double(duration / timeStep) + 0.5;
        }
    };
    /**
     * Refer to [1] for a general review
     * Original literature of RPMD is [3]
     * 
     * Reference:
     * [1] Annual Review of Physical Chemistry 64, 387-413 (2013); https://doi.org/10.1146/annurev-physchem-040412-110122
     * [2] Jos Thijssen. Computational Physics[M].London: Cambridge university press, 2013:197-211
     * [3] J. Chem. Phys. 121, 3368 (2004); https://doi.org/10.1063/1.1777575
     */
    template<class ScalarType,
             unsigned int Dim = 3,
             size_t NumReplica = Dynamic,
             class ForceMatrixAllocator = Utils::HostAllocator<ScalarType>>
    class RPMD final : public RPMDBase<ScalarType> {
        using This = RPMD<ScalarType, Dim, NumReplica, ForceMatrixAllocator>;
        using Base = RPMDBase<ScalarType>;
        using typename Base::PlainScalar;
    public:
        using RingPolymerType = RingPolymer<ScalarType, Dim, NumReplica>;
        using PhaseMatrix = typename RingPolymerType::PhaseMatrix;
        using ForceMatrix = DenseMatrix<ScalarType, MatrixOption::Column | MatrixOption::Vector, Dynamic, Dynamic, ForceMatrixAllocator>;
        using MDCellType = typename RingPolymerType::MDCellType;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using PositionMatrix = typename MDCellType::PositionMatrix;
    private:
        using FireModelType = FireModel<ScalarType, Dim>;

        MDCellType cell;
        RingPolymerType ringPolymer;
        ForceMatrix forceBuffer;

        FFT<ScalarType, 1> fftContract;
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
        template<class ForceModel, class Executor> void updateForce(ForceModel& model);
        template<class KineticModel,
                 class ForceModel,
                 class Executor>
        void nve_step(KineticModel& kineticModel, ForceModel& forceModel);
        template<class KineticModel,
                 class ForceModel,
                 class Executor>
        void nve_step_for(ScalarType duration, KineticModel& kineticModel, ForceModel& forceModel);

        template<class Thermostat,
                 class RandomPoolType,
                 class KineticModel,
                 class ForceModel,
                 class Executor>
        void nvt_step(const Thermostat& thermostat, RandomPoolType& pool, KineticModel& kineticModel, ForceModel& forceModel);
        template<class Thermostat,
                 class RandomPoolType,
                 class KineticModel,
                 class ForceModel,
                 class Executor>
        void nvt_step_for(ScalarType duration, const Thermostat& thermostat, RandomPoolType& pool, KineticModel& kineticModel, ForceModel& forceModel);

        template<class Thermostat,
                 class RandomPoolType,
                 class Barostat,
                 class KineticModel,
                 class ForceModel,
                 class Executor>
        void npt_step(const Thermostat& thermostat, RandomPoolType& pool, Barostat& barostat, KineticModel& kineticModel, ForceModel& forceModel);
        template<class Thermostat,
                 class RandomPoolType,
                 class Barostat,
                 class KineticModel,
                 class ForceModel,
                 class Executor>
        void npt_step_for(ScalarType duration, const Thermostat& thermostat, RandomPoolType& pool, Barostat& barostat, KineticModel& kineticModel, ForceModel& forceModel);

        template<class KineticModel, class ForceModel, class Executor>
        void fire_vstep(FireModelType& fire, KineticModel& kineticModel, ForceModel& forceModel);
        template<BaroType Type, class KineticModel, class ForceModel, class Executor>
        void fire_pstep(CFireModel<ScalarType, Dim, Type>& cfire, KineticModel& kineticModel, ForceModel& forceModel);

        template<class KineticModel, class RandomGenerator> inline void initMomentum(RandomGenerator& gen);
        template<class KineticModel> inline void scaleVelocity();
        void normalizeCentroid();
        [[nodiscard]] MDCellType phaseToCell(size_t replica) const;
        [[nodiscard]] MDCellType contractToCell(size_t contract) const;
        [[nodiscard]] MDCellType makeAverageCell() const;

        template<class KineticModel> [[nodiscard]] ScalarType calcKinetic() const;
        template<class KineticModel> [[nodiscard]] ScalarType calcKinetic(size_t dofIndex) const;
        template<class KineticModel> [[nodiscard]] ScalarType calcKineticPrim() const;
        template<class KineticModel> [[nodiscard]] ScalarType calcKineticPrim(size_t dofIndex) const;
        [[nodiscard]] inline ScalarType calcKineticClassical() const;

        template<class ForceModel, class Executor> [[nodiscard]] ScalarType calcPotential(const ForceModel& model) const;
        template<class ForceModel> [[nodiscard]] ScalarType calcPotentialClassical(const ForceModel& model) const;

        [[nodiscard]] ScalarType calcClassicalElastic() const;
        template<class ForceModel> [[nodiscard]] ScalarType calcClassicalInternalEnergy(const ForceModel& model) const;

        template<class KineticModel> [[nodiscard]] ScalarType calcTemperature() const;
        template<class KineticModel, class ForceModel, class Executor>
        [[nodiscard]] ScalarType calcPressThermo(ForceModel& model) const;

        template<class ForceModel, class Executor>
        [[nodiscard]] LatticeMatrix makeStressPrim(ForceModel& model) const;
        template<class ForceModel, class Executor>
        [[nodiscard]] LatticeMatrix makeStressVirial(ForceModel& model) const;
        template<class ForceModel, class Executor>
        [[nodiscard]] LatticeMatrix makeStressClassical(ForceModel& model) const;

        template<class KineticModel, class ForceModel, class Executor>
        Vector<ScalarType> testNVE(ScalarType duration, KineticModel& kineticModel, ForceModel& forceModel) const;

        void read(const H5Location& loc, const char* name);
        void write(H5Location& loc, const char* name) const;
        void swap(RPMD& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] constexpr unsigned int getDim() const noexcept { return Dim; }
        [[nodiscard]] const LatticeMatrix& getLattice() const noexcept { return cell.getLattice(); }
        [[nodiscard]] const typename MDCellType::InvLatticeMatrix& getInvLattice() const noexcept { return cell.getInvLattice(); }
        [[nodiscard]] const typename MDCellType::MassVector& getMassVec() const noexcept { return cell.getMassVec(); }
        [[nodiscard]] size_t getNumParticle() const noexcept { return cell.getNumParticle(); }
        [[nodiscard]] ScalarType getVolume() const noexcept { return cell.getVolume(); }
        [[nodiscard]] RingPolymerType& getRingPolymer() noexcept { return ringPolymer; }
        [[nodiscard]] const RingPolymerType& getRingPolymer() const noexcept { return ringPolymer; }
        [[nodiscard]] PhaseMatrix& getPhaseMatrix() noexcept { return ringPolymer.asMatrix(); }
        [[nodiscard]] const PhaseMatrix& getPhaseMatrix() const noexcept { return ringPolymer.asMatrix(); }
        [[nodiscard]] size_t getDOF() const noexcept { return ringPolymer.getDOF(); }
        [[nodiscard]] size_t getNumReplica() const noexcept { return ringPolymer.getNumReplica(); }
        [[nodiscard]] size_t getKSpaceSize() const noexcept { return ringPolymer.getKSpaceSize(); }
        [[nodiscard]] const ForceMatrix& getForce() const noexcept { return forceBuffer; }
        [[nodiscard]] size_t getNumContract() const noexcept { return fftContract.getRSpaceSize(); }
        [[nodiscard]] bool isContractEnabled() const noexcept { return (NumReplica != 1) && (getNumReplica() != getNumContract()); }
        [[nodiscard]] ScalarType getTemperature() const noexcept { return temperatureT; }
        [[nodiscard]] ScalarType getTimeStep() const noexcept { return timeStep; }
        /* Setters */
        void setLattice(LatticeMatrix lattice) { cell.setLattice(std::move(lattice)); }
        inline void setTemperature(ScalarType temperature);
        void setTimeStep(ScalarType timeStep_);
    private:
        void toContractBeadRepr(size_t posID);
        void forceToNormRepr(size_t posID);
        void forceToBeadRepr(size_t posID);
        void contract();
        void decontract();
        void forceStep(ScalarType deltaT);
        bool checkCentroid() const;
        void checkParam() const;
    };
}

#include "MDImpl/RPMDImpl.h"
