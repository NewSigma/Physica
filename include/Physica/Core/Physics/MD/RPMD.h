/*
 * Copyright 2022 WeiBo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Flatten.h"
#include "Physica/Core/Math/Statistics/NumCharacter.h"
#include "Physica/Core/Physics/PhyConst.h"
#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"
#include "MDCell.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] M. Ceriotti, M. Parrinello, T. E. Markland and D. E. Manolopoulos, J. Chem. Phys. 133, 124104 (2010).
     * [2] Habershon S, Manolopoulos D E, Markland T E, et al. Ring-Polymer Molecular Dynamics: Quantum Effects in Chemical Dynamics from Classical Trajectories in an Extended Phase Space[J]. Annual Review of Physical Chemistry, 2013, 64(1):387-413.
     * [3] Rossi M, Ceriotti M, Manolopoulos D E. How to remove the spurious resonances from ring polymer molecular dynamics[J]. J. Chem. Phys, 2014, 140(23):5106.
     * [4] Jos Thijssen. Computational Physics[M].London: Cambridge university press, 2013:197-211
     * [5] Liu J, Li D, Liu X. A simple and accurate algorithm for path integral molecular dynamics with the Langevin thermostat[J]. J. Chem. Phys, 2016, 145(2):1291-1301.
     * [6] T. E. Markland, D. E. Manolopoulos. An efficient ring polymer contraction scheme for imaginary time path integral simulations[J]. J. Chem. Phys. 129, 024105 (2008)
     * 
     * TODO: replace several ScalarType to PosScalarType
     */
    template<class ScalarType, class PosScalarType>
    class RPMD final {
        using BufferType = DenseMatrix<ComplexScalar<ScalarType>, MatrixOption::Row | MatrixOption::Vector, 2>;
    public:
        using PhaseMatrixType = DenseMatrix<PosScalarType, MatrixOption::Row | MatrixOption::Vector>;
        using ForceMatrix = DenseMatrix<ScalarType, MatrixOption::Column | MatrixOption::Vector>;
        using MDCellType = MDCell<ScalarType, PosScalarType>;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using PositionMatrix = typename MDCellType::PositionMatrix;
        constexpr static unsigned int Dim = MDCellType::Dim;
    private:
        MDCellType cell;
        FFT<ScalarType, 1> fft;
        PhaseMatrixType phaseMatrix;
        ForceMatrix forceBuffer;

        FFT<PosScalarType, 1> fftContract;
        PhaseMatrixType posContract;
        ForceMatrix forceContract;

        BufferType buffer;
        /* Constant */
        ScalarType temperatureT;
        ScalarType thermostatTime;
        ScalarType timeStep;
        ScalarType repBeta;
        ScalarType omegaW;
    public:
        RPMD(MDCellType cell_,
             size_t numReplica,
             size_t numContract,
             ScalarType temperatureT_,
             ScalarType thermostatTime_,
             ScalarType timeStep_);
        RPMD(const RPMD&) = default;
        RPMD(RPMD&&) noexcept = default;
        ~RPMD() = default;
        /* Operations */
        template<class ForceModel, class Executor>
        void updateForce(const ForceModel& model);
        template<class RandomGenerator, class ForceModel, class Executor>
        void nvt_step(RandomGenerator& gen, const ForceModel& model);
        template<class ForceModel, class Executor>
        void nve_step(const ForceModel& model);
        template<class RandomGenerator, class Barostat, class ForceModel, class Executor>
        void npt_step(RandomGenerator& gen, Barostat& barostat, const ForceModel& model);
        template<class RandomGenerator, class ForceModel, class Executor>
        void nvt_step_for(ScalarType duration, RandomGenerator& gen, const ForceModel& model);
        template<class ForceModel, class Executor>
        void nve_step_for(ScalarType duration, const ForceModel& model);
        template<class RandomGenerator, class Barostat, class ForceModel, class Executor>
        void npt_step_for(ScalarType duration, RandomGenerator& gen, Barostat& barostat, const ForceModel& model);
        template<class RandomGenerator>
        void initMomentum(RandomGenerator& gen);
        void removeDrift();
        void scaleVelocity();
        void normalizeCentroid();
        [[nodiscard]] MDCellType phaseToCell(size_t replica) const;
        [[nodiscard]] MDCellType contractToCell(size_t contract) const;
        void checkParam() const;
        /* Getters */
        [[nodiscard]] constexpr unsigned int getDim() const noexcept { return Dim; }
        [[nodiscard]] const LatticeMatrix& getLattice() const noexcept { return cell.getLattice(); }
        [[nodiscard]] const typename MDCellType::InvLatticeMatrix& getInvLattice() const noexcept { return cell.getInvLattice(); }
        [[nodiscard]] const typename MDCellType::MassVector& getMassVec() const noexcept { return cell.getMassVec(); }
        [[nodiscard]] size_t getNumParticle() const noexcept { return cell.getNumParticle(); }
        [[nodiscard]] PosScalarType getVolume() const noexcept { return cell.getVolume(); }
        [[nodiscard]] const PhaseMatrixType& getPhaseMatrix() const noexcept { return phaseMatrix; }
        [[nodiscard]] PhaseMatrixType& getPhaseMatrix() noexcept { return phaseMatrix; }
        [[nodiscard]] size_t getNumReplica() const noexcept { return phaseMatrix.getColumn(); }
        [[nodiscard]] size_t getNumContract() const noexcept { return fftContract.getRSpaceSize(); }
        [[nodiscard]] bool isContractEnabled() const noexcept { return getNumReplica() != getNumContract(); }
        [[nodiscard]] size_t getDOF() const noexcept { return Dim * getNumParticle(); }
        [[nodiscard]] ScalarType getTemperature() const noexcept { return temperatureT; }
        [[nodiscard]] ScalarType getOmegaW() const noexcept { return omegaW; }

        [[nodiscard]] PositionMatrix makeCentroidPos() const;
        [[nodiscard]] MDCellType makeAverageCell() const;
        [[nodiscard]] PositionMatrix makeCentroidMomentum() const;
        [[nodiscard]] PositionMatrix getMomentum(size_t replica) const;
        [[nodiscard]] const ForceMatrix& getForce() const noexcept { return forceBuffer; }

        [[nodiscard]] ScalarType getClassicalKinetic() const;
        template<class ForceModel>
        [[nodiscard]] ScalarType getClassicalPotentialEnergy(const ForceModel& model) const;
        [[nodiscard]] ScalarType getClassicalElastic() const;
        template<class ForceModel>
        [[nodiscard]] ScalarType getClassicalInternalEnergy(const ForceModel& model) const;
        [[nodiscard]] ScalarType calcKinetic() const;
        template<class ForceModel, class Executor>
        [[nodiscard]] ScalarType calcPotential(const ForceModel& model) const;
        [[nodiscard]] ScalarType calcTemperature() const;
        [[nodiscard]] LatticeMatrix makeStress() const;
        /* Setters */
        void setTemperature(ScalarType temperature);
        void setThermostatTime(ScalarType time) { thermostatTime = time; }
    private:
        void toNormalRepr(size_t posID);
        void toBeadRepr(size_t posID);
        void toContractBeadRepr(size_t posID);
        void forceToNormRepr(size_t posID);
        void forceToBeadRepr(size_t posID);
        void contract();
        void decontract();
        template<class RandomGenerator>
        void thermostatStep(RandomGenerator& gen, ScalarType deltaT);
        void thermostatImpl(size_t mode_index, ScalarType deltaT, ScalarType viscosityY, ScalarType factor, ComplexScalar<ScalarType> random);
        void forceStep(ScalarType deltaT);
        void dynamicStep(ScalarType deltaT);
        template<class Barostat>
        void npt_dynamicStep(Barostat& barostat, ScalarType deltaT);
        bool checkCentroid() const;
    };
}

#include "MDImpl/RPMDImpl.h"
