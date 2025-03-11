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

#include "Physica/Core/Math/Transform/FFT.h"
#include "MDImpl/RingPolymer.h"
#include "Barostat/BaroType.h"

namespace Physica {
    template<Scalar T, unsigned int Dim> class FireModel;
    template<Scalar T, unsigned int Dim, BaroType Type> class CFireModel;
    /**
     * Refer to [1] for a general review
     * Original literature of RPMD is [3]
     * 
     * Reference:
     * [1] Annual Review of Physical Chemistry 64, 387-413 (2013); https://doi.org/10.1146/annurev-physchem-040412-110122
     * [2] J. H. Thijssen. Computational Physics[M]. London: Cambridge University Press, 2013:197-211
     * [3] J. Chem. Phys. 121, 3368 (2004); https://doi.org/10.1063/1.1777575
     */
    template<Scalar T,
             unsigned int Dim = 3,
             size_t NumReplica = Dynamic,
             class ForceMatrixAllocator = HostAllocator<T>>
    class RPMD final {
        using This = RPMD<T, Dim, NumReplica, ForceMatrixAllocator>;
        using Tv = T::ValueType;
    public:
        using RingPolymerType = RingPolymer<T, Dim, NumReplica>;
        using PhaseMatrix = RingPolymerType::PhaseMatrix;
        using ForceMatrix = DenseMatrix<T, MatrixOption::Col | MatrixOption::Vector, Dynamic, Dynamic, ForceMatrixAllocator>;
        using MDCellType = RingPolymerType::MDCellType;
        using LatticeMatrix = MDCellType::LatticeMatrix;
        using PositionMatrix = MDCellType::PositionMatrix;
    private:
        using FireModelType = FireModel<T, Dim>;

        MDCellType cell;
        RingPolymerType ringPolymer;
        ForceMatrix forceBuffer;

        FFT<T, 1> fftContract;
        PhaseMatrix posContract;
        ForceMatrix forceContract;
        /* Constant */
        T temperatureT;
        T timeStep;
    public:
        RPMD(MDCellType cell_,
             size_t numReplica,
             size_t numContract,
             T temperatureT_,
             T timeStep_);
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
        void nve_step_for(T duration, KineticModel& kineticModel, ForceModel& forceModel);

        template<class Thermostat,
                 RNG R,
                 class KineticModel,
                 class ForceModel,
                 class Executor>
        void nvt_step(const Thermostat& thermostat, KineticModel& kineticModel, ForceModel& forceModel);
        template<class Thermostat,
                 RNG R,
                 class KineticModel,
                 class ForceModel,
                 class Executor>
        void nvt_step_for(T duration, const Thermostat& thermostat, KineticModel& kineticModel, ForceModel& forceModel);

        template<class Thermostat,
                 RNG R,
                 class Barostat,
                 class KineticModel,
                 class ForceModel,
                 class Executor>
        void npt_step(const Thermostat& thermostat, Barostat& barostat, KineticModel& kineticModel, ForceModel& forceModel);
        template<class Thermostat,
                 RNG R,
                 class Barostat,
                 class KineticModel,
                 class ForceModel,
                 class Executor>
        void npt_step_for(T duration, const Thermostat& thermostat, Barostat& barostat, KineticModel& kineticModel, ForceModel& forceModel);

        template<class KineticModel, class ForceModel, class Executor>
        void fire_vstep(FireModelType& fire, KineticModel& kineticModel, ForceModel& forceModel);
        template<BaroType Type, class KineticModel, class ForceModel, class Executor>
        void fire_pstep(CFireModel<T, Dim, Type>& cfire, KineticModel& kineticModel, ForceModel& forceModel);

        template<class KineticModel, RNG R> inline void initMomentum();
        template<class KineticModel> inline void scaleVelocity();
        void normalizeCentroid();
        [[nodiscard]] MDCellType phaseToCell(size_t replica) const;
        [[nodiscard]] MDCellType contractToCell(size_t contract) const;
        [[nodiscard]] MDCellType makeAverageCell() const;

        template<class KineticModel> [[nodiscard]] T calcKinetic() const;
        template<class KineticModel> [[nodiscard]] T calcKinetic(size_t dofIndex) const;
        template<class KineticModel> [[nodiscard]] T calcKineticPrim() const;
        template<class KineticModel> [[nodiscard]] T calcKineticPrim(size_t dofIndex) const;
        [[nodiscard]] inline T calcKineticClassical() const;

        template<class ForceModel, class Executor> [[nodiscard]] T calcPotential(const ForceModel& model) const;
        template<class ForceModel> [[nodiscard]] T calcPotentialClassical(const ForceModel& model) const;

        [[nodiscard]] T calcClassicalElastic() const;
        template<class ForceModel> [[nodiscard]] T calcClassicalInternalEnergy(const ForceModel& model) const;

        template<class KineticModel> [[nodiscard]] T calcTemperature() const;
        template<class KineticModel, class ForceModel, class Executor>
        [[nodiscard]] T calcPressThermo(ForceModel& model) const;

        template<class ForceModel, class Executor>
        [[nodiscard]] LatticeMatrix makeStressPrim(ForceModel& model) const;
        template<class ForceModel, class Executor>
        [[nodiscard]] LatticeMatrix makeStressVirial(ForceModel& model) const;
        template<class ForceModel, class Executor>
        [[nodiscard]] LatticeMatrix makeStressClassical(ForceModel& model) const;

        template<class KineticModel, class ForceModel, class Executor>
        VectorND<T> testNVE(T duration, KineticModel& kineticModel, ForceModel& forceModel) const;

        void read(const H5Loc& loc, const char* name);
        void write(H5Loc& loc, const char* name) const;
        void swap(RPMD& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] constexpr unsigned int getDim() const noexcept { return Dim; }
        [[nodiscard]] const LatticeMatrix& getLattice() const noexcept { return cell.getLattice(); }
        [[nodiscard]] const MDCellType::InvLatticeMatrix& getInvLattice() const noexcept { return cell.getInvLattice(); }
        [[nodiscard]] const MDCellType::MassVector& getMassVec() const noexcept { return cell.getMassVec(); }
        [[nodiscard]] size_t getNumParticle() const noexcept { return cell.getNumParticle(); }
        [[nodiscard]] T getVolume() const noexcept { return cell.getVolume(); }
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
        [[nodiscard]] T getTemperature() const noexcept { return temperatureT; }
        [[nodiscard]] T getTimeStep() const noexcept { return timeStep; }
        /* Setters */
        void setLattice(LatticeMatrix lattice) { cell.setLattice(std::move(lattice)); }
        inline void setTemperature(T temperature);
        void setTimeStep(T timeStep_);
        /* Static members */
        [[nodiscard]] static uint64_t durationToStep(Tv duration, Tv timeStep);
    private:
        void toContractBeadRepr(size_t posID);
        void forceToNormRepr(size_t posID);
        void forceToBeadRepr(size_t posID);
        void contract();
        void decontract();
        void forceStep(T deltaT);
        bool checkCentroid() const;
        void checkParam() const;
    };
}

#include "MDImpl/RPMDImpl.h"
