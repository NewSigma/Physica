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
        using ForceMatrix = DenseMatrix<T, MatrixOption::Col, Dynamic, Dynamic, ForceMatrixAllocator>;
        using MDCellType = RingPolymerType::MDCellType;
        using LatticeMatrix = MDCellType::LatticeMatrix;
        using PositionMatrix = MDCellType::PositionMatrix;
    private:
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
        RPMD() = default;
        RPMD(MDCellType cell_,
             size_t numReplica,
             size_t numContract,
             T temperatureT_,
             T timeStep_);
        RPMD(const RPMD&) = default;
        RPMD(RPMD&&) noexcept = default;
        ~RPMD() = default;
        /* Operators */
        RPMD& operator=(RPMD obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<ExecutePolicy P>
        void updateForce(auto& model);

        template<ExecutePolicy P = Sequential>
        void nve_step(auto& kineticModel, auto& forceModel);
        template<ExecutePolicy P = Sequential>
        void nve_step_for(T duration, auto& kineticModel, auto& forceModel);

        template<RNG R, ExecutePolicy P = Sequential>
        void nvt_step(const auto& thermostat, auto& kineticModel, auto& forceModel);
        template<RNG R, ExecutePolicy P = Sequential>
        void nvt_step_for(T duration, const auto& thermostat, auto& kineticModel, auto& forceModel);

        template<RNG R, ExecutePolicy P = Sequential>
        void npt_step(const auto& thermostat, auto& barostat, auto& kineticModel, auto& forceModel);
        template<RNG R, ExecutePolicy P = Sequential>
        void npt_step_for(T duration, const auto& thermostat, auto& barostat, auto& kineticModel, auto& forceModel);

        template<ExecutePolicy P = Sequential>
        void fire_vstep(FireModel<T, Dim>& fire, auto& kineticModel, auto& forceModel);
        template<BaroType Type, ExecutePolicy P = Sequential>
        void fire_pstep(CFireModel<T, Dim, Type>& cfire, auto& kineticModel, auto& forceModel);

        template<class KineticModel, RNG R>
        void initMomentum();
        template<class KineticModel>
        void scaleVelocity();
        void normalizeCentroid();
        [[nodiscard]] MDCellType phaseToCell(size_t replica) const;
        [[nodiscard]] MDCellType contractToCell(size_t contract) const;
        [[nodiscard]] MDCellType makeAverageCell() const;

        template<class KineticModel> [[nodiscard]] T calcKinetic() const;
        template<class KineticModel> [[nodiscard]] T calcKinetic(size_t dofIndex) const;
        template<class KineticModel> [[nodiscard]] T calcKineticPrim() const;
        template<class KineticModel> [[nodiscard]] T calcKineticPrim(size_t dofIndex) const;
        [[nodiscard]] T calcKineticClassical() const;

        template<ExecutePolicy P = Sequential>
        [[nodiscard]] T calcPotential(const auto& model) const;
        [[nodiscard]] T calcPotentialClassical(const auto& model) const;

        [[nodiscard]] T calcClassicalElastic() const;
        [[nodiscard]] T calcClassicalInternalEnergy(const auto& forceModel) const;

        template<class KineticModel> [[nodiscard]] T calcTemperature() const;
        template<class KineticModel, ExecutePolicy P>
        [[nodiscard]] T calcPressThermo(auto& model) const;

        template<ExecutePolicy P>
        [[nodiscard]] LatticeMatrix makeStressPrim(auto& forceModel) const;
        template<ExecutePolicy P>
        [[nodiscard]] LatticeMatrix makeStressVirial(auto& forceModel) const;
        template<ExecutePolicy P>
        [[nodiscard]] LatticeMatrix makeStressClassical(auto& forceModel) const;

        template<ExecutePolicy P = Sequential>
        [[nodiscard]] VectorND<T> testNVE(T duration, auto& kineticModel, auto& forceModel) const;

        void read(const H5Loc& loc, const char* name);
        void write(H5Loc& loc, const char* name) const;
        void swap(RPMD& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] constexpr unsigned int getDim() const noexcept { return Dim; }
        [[nodiscard]] const auto& getLattice() const noexcept { return cell.getLattice(); }
        [[nodiscard]] const auto& getInvLattice() const noexcept { return cell.getInvLattice(); }
        [[nodiscard]] const auto& getMassVec() const noexcept { return cell.getMassVec(); }
        [[nodiscard]] size_t getNumParticle() const noexcept { return cell.getNumParticle(); }
        [[nodiscard]] T getVolume() const noexcept { return cell.getVolume(); }
        [[nodiscard]] auto&& getRingPolymer(this auto&&) noexcept;
        [[nodiscard]] auto&& getPhaseMatrix(this auto&&) noexcept;
        [[nodiscard]] size_t getDOF() const noexcept { return ringPolymer.getDOF(); }
        [[nodiscard]] size_t getNumReplica() const noexcept { return ringPolymer.getNumReplica(); }
        [[nodiscard]] size_t getKSpaceSize() const noexcept { return ringPolymer.getKSpaceSize(); }
        [[nodiscard]] const auto& getForce() const noexcept { return forceBuffer; }
        [[nodiscard]] size_t getNumContract() const noexcept { return fftContract.getRSpaceSize(); }
        [[nodiscard]] bool isContractEnabled() const noexcept { return (NumReplica != 1) && (getNumReplica() != getNumContract()); }
        [[nodiscard]] T getTemperature() const noexcept { return temperatureT; }
        [[nodiscard]] T getTimeStep() const noexcept { return timeStep; }
        /* Setters */
        void setLattice(LatticeMatrix lattice) { cell.setLattice(std::move(lattice)); }
        void setTemperature(T temperature) noexcept;
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
        template<class KineticModel, class ForceModel>
        constexpr static void checkRelaxParam();
    };
}

#include "MDImpl/RPMDImpl.h"
