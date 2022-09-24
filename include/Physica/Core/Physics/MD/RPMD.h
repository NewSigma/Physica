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
     * [3] Rossi M, Ceriotti M, Manolopoulos D E. How to remove the spurious resonances from ring polymer molecular dynamics[J]. Journal of Chemical Physics, 2014, 140(23):5106.
     * [4] Jos Thijssen. Computational Physics[M].London: Cambridge university press, 2013:197-211
     * [5] Liu J, Li D, Liu X. A simple and accurate algorithm for path integral molecular dynamics with the Langevin thermostat[J]. The Journal of Chemical Physics, 2016, 145(2):1291-1301.
     */
    template<class ScalarType, class PosScalarType>
    class RPMD final {
        using BufferType = DenseMatrix<ComplexScalar<ScalarType>, MatrixOption::Row | MatrixOption::Vector, 2>;
    public:
        using PhasePosType = DenseMatrix<PosScalarType, MatrixOption::Row | MatrixOption::Vector>;
        using ForceMatrix = DenseMatrix<ScalarType, MatrixOption::Column | MatrixOption::Vector>;
        using MDCellType = MDCell<ScalarType, PosScalarType>;
        using PositionMatrix = typename MDCellType::PositionMatrix;
        constexpr static unsigned int Dim = MDCellType::Dim;
    private:
        MDCellType cell;
        FFT<ScalarType, 1> fft;
        PhasePosType phasePosX;
        ForceMatrix forceBuffer;
        BufferType buffer;
        ScalarType temperatureT;
        ScalarType thermostatTime;
        ScalarType timeStep;

        ScalarType repBeta;
        ScalarType omegaW;
    public:
        RPMD(MDCellType cell_, size_t numReplica, ScalarType temperatureT_, ScalarType thermostatTime_, ScalarType timeStep_);
        /* Operations */
        template<class ForceModel, class Executor = Parallel::SequentialExecutor>
        void updateForce(const ForceModel& model);
        template<class RandomGenerator, class ForceModel, class Executor = Parallel::SequentialExecutor>
        void nvt_step(RandomGenerator& gen, const ForceModel& force);
        template<class ForceModel, class Executor = Parallel::SequentialExecutor>
        void nve_step(const ForceModel& force);
        template<class RandomGenerator, class ForceModel, class Executor = Parallel::SequentialExecutor>
        void nvt_step_for(ScalarType duration, RandomGenerator& gen, const ForceModel& force);
        template<class ForceModel, class Executor = Parallel::SequentialExecutor>
        void nve_step_for(ScalarType duration, const ForceModel& force);
        template<class RandomGenerator>
        void initMomentum(RandomGenerator& gen);
        void removeDrift();
        void scaleVelocity();
        template<class RandomGenerator, class ForceModel, class Executor = Parallel::SequentialExecutor>
        [[nodiscard]] bool isStableNVT(size_t numStep, RandomGenerator& gen, const ForceModel& force, double precision);
        /* Getters */
        [[nodiscard]] const typename MDCellType::LatticeMatrix& getLattice() const noexcept { return cell.getLattice(); }
        [[nodiscard]] const typename MDCellType::MassVector& getMassVec() const noexcept { return cell.getMassVec(); }
        [[nodiscard]] size_t getNumParticle() const noexcept { return cell.getNumParticle(); }
        [[nodiscard]] const PhasePosType& getPhasePos() const noexcept { return phasePosX; }
        [[nodiscard]] PhasePosType& getPhasePos() noexcept { return phasePosX; }
        [[nodiscard]] size_t getNumReplica() const noexcept { return phasePosX.getColumn(); }
        [[nodiscard]] MDCellType phaseToCell(size_t replica) const;
        [[nodiscard]] size_t getDOF() const noexcept { return Dim * getNumParticle(); }
        [[nodiscard]] ScalarType getTemperature() const noexcept { return temperatureT; }
        [[nodiscard]] PeriodicCell<PosScalarType, 3> makeAverageCell() const;
        [[nodiscard]] PositionMatrix getPos() const;
        [[nodiscard]] PositionMatrix getMomentum() const;
        [[nodiscard]] const ForceMatrix& getForce() const noexcept { return forceBuffer; }

        [[nodiscard]] ScalarType getClassicalKinetic() const;
        template<class ForceModel>
        [[nodiscard]] ScalarType getClassicalPotentialEnergy(ForceModel model) const;
        [[nodiscard]] ScalarType getClassicalElastic() const;
        template<class ForceModel>
        [[nodiscard]] ScalarType getClassicalInternalEnergy(ForceModel model) const;
        template<class ForceModel>
        [[nodiscard]] ScalarType computeKinetic(ForceModel model) const;
        [[nodiscard]] ScalarType estimateTemperature() const;
        /* Setters */
        void setTemperature(ScalarType temperature);
        void setThermostatTime(ScalarType time) { thermostatTime = time; }
    private:
        void toNormalRepr(size_t posID);
        void toBeadRepr(size_t posID);
        template<class RandomGenerator>
        void thermostatStep(RandomGenerator& gen, ScalarType deltaT);
        void thermostatImpl(size_t mode_index, ScalarType deltaT, ScalarType viscosityY, ScalarType factor, ComplexScalar<ScalarType> random);
        void forceStep(ScalarType deltaT);
        void dynamicStep(ScalarType deltaT);
        void normalizeCentroid();
        bool checkCentroid() const;
    };

    template<class ScalarType, class PosScalarType>
    RPMD<ScalarType, PosScalarType>::RPMD(MDCellType cell_, size_t numReplica, ScalarType temperatureT_, ScalarType thermostatTime_, ScalarType timeStep_)
            : cell(std::move(cell_))
            , fft(numReplica, 1)
            , thermostatTime(std::move(thermostatTime_))
            , timeStep(std::move(timeStep_)) {
        const size_t dof = getDOF();
        phasePosX.resize(2 * dof, numReplica);
        forceBuffer.resize(dof, numReplica);
        buffer.resize(2, fft.getKSpaceSize());

        auto momentum = phasePosX.topRows(dof);
        momentum = ScalarType::Zero();

        /* Fill pos */ {
            size_t index = dof;
            for (auto elem : cell.getPos()) {
                phasePosX(index, 0) = elem;
                ++index;
            }
            for (size_t i = 1; i < getNumReplica(); ++i) {
                auto phasePos = phasePosX.col(i);
                auto pos = phasePos.tail(dof);
                pos = phasePosX.col(0).tail(dof);
            }
        }
        setTemperature(temperatureT_);
    }

    template<class ScalarType, class PosScalarType>
    template<class ForceModel, class Executor>
    void RPMD<ScalarType, PosScalarType>::updateForce(const ForceModel& model) {
        auto kernel = [&](unsigned int replica) {
            MDCellType cell = phaseToCell(replica);
            cell.normalizeCell();
            auto saveTo = forceBuffer.col(replica);
            saveTo = model.force(std::move(cell));
        };
        Executor::parallel_for(kernel, getNumReplica(), Executor::getNumThread()).wait();
    }

    template<class ScalarType, class PosScalarType>
    template<class RandomGenerator, class ForceModel, class Executor>
    void RPMD<ScalarType, PosScalarType>::nvt_step(RandomGenerator& gen, const ForceModel& force) {
        forceStep(timeStep * 0.5);
        dynamicStep(timeStep * 0.5);
        thermostatStep(gen, timeStep);
        dynamicStep(timeStep * 0.5);
        updateForce<ForceModel, Executor>(force);
        forceStep(timeStep * 0.5);
    }

    template<class ScalarType, class PosScalarType>
    template<class ForceModel, class Executor>
    void RPMD<ScalarType, PosScalarType>::nve_step(const ForceModel& force) {
        forceStep(timeStep * 0.5);
        dynamicStep(timeStep);
        normalizeCentroid();
        updateForce<ForceModel, Executor>(force);
        forceStep(timeStep * 0.5);
    }

    template<class ScalarType, class PosScalarType>
    template<class RandomGenerator, class ForceModel, class Executor>
    void RPMD<ScalarType, PosScalarType>::nvt_step_for(ScalarType duration, RandomGenerator& gen, const ForceModel& force) {
        uint64_t step = double(duration / timeStep) + 0.5;
        for (uint64_t _ = 0; _ < step; ++_)
            nvt_step<RandomGenerator, ForceModel, Executor>(gen, force);
    }

    template<class ScalarType, class PosScalarType>
    template<class ForceModel, class Executor>
    void RPMD<ScalarType, PosScalarType>::nve_step_for(ScalarType duration, const ForceModel& force) {
        uint64_t step = double(duration / timeStep) + 0.5;
        for (uint64_t _ = 0; _ < step; ++_)
            nve_step<ForceModel, Executor>(force);
    }

    template<class ScalarType, class PosScalarType>
    template<class RandomGenerator>
    void RPMD<ScalarType, PosScalarType>::initMomentum(RandomGenerator& gen) {
        std::normal_distribution<> dist{};
        const size_t dof = getDOF();
        Vector<ScalarType, 3> driftMomentum(3, 0);
        for (size_t i = 0; i < dof; ++i) {
            const auto mass = cell.getMass(i / Dim);
            const size_t direction = i % Dim;
            const ScalarType factor = sqrt(repBeta * mass);
            for (size_t j = 0; j < getNumReplica(); ++j) {
                const ScalarType temp = factor * dist(gen);
                phasePosX(i, j) = temp;
                driftMomentum[direction] += temp;
            }
        }
        driftMomentum *= Core::reciprocal(ScalarType(getNumParticle() * getNumReplica()));

        for (size_t i = 0; i < dof; ++i) {
            auto row = phasePosX.row(i);
            row -= driftMomentum[i % 3];
        }
    }

    template<class ScalarType, class PosScalarType>
    void RPMD<ScalarType, PosScalarType>::removeDrift() {
        const size_t dof = getDOF();
        Vector<ScalarType, 3> driftMomentum(3, 0);
        for (size_t i = 0; i < dof; ++i) {
            const size_t direction = i % Dim;
            for (size_t j = 0; j < getNumReplica(); ++j)
                driftMomentum[direction] += phasePosX(i, j);
        }
        driftMomentum *= Core::reciprocal(ScalarType(getNumParticle() * getNumReplica()));

        for (size_t i = 0; i < dof; ++i) {
            auto row = phasePosX.row(i);
            row -= driftMomentum[i % 3];
        }
    }

    template<class ScalarType, class PosScalarType>
    void RPMD<ScalarType, PosScalarType>::scaleVelocity() {
        const ScalarType temperatureNow = estimateTemperature();
        assert(temperatureNow.isPositive());
        const size_t dof = getDOF();
        const ScalarType factor = sqrt(temperatureT / temperatureNow);
        auto momentum = phasePosX.topRows(dof);
        momentum *= factor;
    }

    template<class ScalarType, class PosScalarType>
    template<class RandomGenerator, class ForceModel, class Executor>
    bool RPMD<ScalarType, PosScalarType>::isStableNVT(size_t numStep, RandomGenerator& gen, const ForceModel& force, double precision) {
        ScalarType kinetic = 0;
        ScalarType squared_kinetic = 0;
        for (size_t i = 0; i < numStep; ++i) {
            nvt_step<RandomGenerator, ForceModel, Executor>(gen, force);
            const ScalarType temp = getClassicalKinetic();
            toNextMean(kinetic, i, temp);
            toNextMean(squared_kinetic, i, square(temp));
        }
        const ScalarType factor = ScalarType(getDOF() * getNumReplica()) / ScalarType(getDOF() * getNumReplica() + 2);
        return scalarNear(square(kinetic), factor * squared_kinetic, precision);
    }

    template<class ScalarType, class PosScalarType>
    typename RPMD<ScalarType, PosScalarType>::MDCellType RPMD<ScalarType, PosScalarType>::phaseToCell(size_t replica) const {
        PositionMatrix pos(getNumParticle(), 3);
        auto phase = phasePosX.col(replica);
        size_t index = getDOF();
        for (auto& elem : pos) {
            elem = PosScalarType(phase[index]);
            ++index;
        }
        return MDCellType(cell.getLattice(), std::move(pos), cell.getMassVec());
    }

    template<class ScalarType, class PosScalarType>
    PeriodicCell<PosScalarType, 3> RPMD<ScalarType, PosScalarType>::makeAverageCell() const {
        return PeriodicCell<PosScalarType, 3>(getLattice(), getPos());
    }

    template<class ScalarType, class PosScalarType>
    typename RPMD<ScalarType, PosScalarType>::PositionMatrix RPMD<ScalarType, PosScalarType>::getPos() const {
        PositionMatrix result(getNumParticle(), 3);
        const size_t dof = getDOF();
        size_t index = dof;
        for (auto& elem : result) {
            elem = PosScalarType(mean(phasePosX.row(index)));
            ++index;
        }
        return result;
    }

    template<class ScalarType, class PosScalarType>
    typename RPMD<ScalarType, PosScalarType>::PositionMatrix RPMD<ScalarType, PosScalarType>::getMomentum() const {
        PositionMatrix result(getNumParticle(), 3, 0);
        size_t index = 0;
        for (auto& elem : result) {
            elem = PosScalarType(mean(phasePosX.row(index)));
            ++index;
        }
        return result;
    }

    template<class ScalarType, class PosScalarType>
    ScalarType RPMD<ScalarType, PosScalarType>::getClassicalKinetic() const {
        const size_t dof = getDOF();
        ScalarType classical_kinetic = 0;
        for (size_t i = 0; i < dof; ++i) {
            const auto mass = cell.getMass(i / Dim);
            auto p = phasePosX.row(i);
            classical_kinetic += square(p.asVector()).sum() / (mass * 2);
        }
        return classical_kinetic;
    }

    template<class ScalarType, class PosScalarType>
    template<class ForceModel>
    ScalarType RPMD<ScalarType, PosScalarType>::getClassicalPotentialEnergy(ForceModel model) const {
        ScalarType result = 0;
        for (size_t i = 0; i < getNumReplica(); ++i)
            result += model.potentialEnergy(phaseToCell(i));
        return result;
    }

    template<class ScalarType, class PosScalarType>
    ScalarType RPMD<ScalarType, PosScalarType>::getClassicalElastic() const {
        const size_t dof = getDOF();
        auto pos = phasePosX.bottomRows(dof);
        ScalarType result = 0;
        for (size_t i = 0; i < dof; ++i) {
            const ScalarType mass = cell.getMass(i / Dim);
            for (size_t j = 0; j < getNumReplica(); ++j)
                result += mass * square(omegaW * (pos(i, j) - pos(i, (j + 1) % getNumReplica()))) * 0.5;
        }
        return result;
    }

    template<class ScalarType, class PosScalarType>
    template<class ForceModel>
    ScalarType RPMD<ScalarType, PosScalarType>::getClassicalInternalEnergy(ForceModel model) const {
        return getClassicalKinetic() + getClassicalPotentialEnergy<ForceModel>(model) + getClassicalElastic();
    }

    template<class ScalarType, class PosScalarType>
    template<class ForceModel>
    ScalarType RPMD<ScalarType, PosScalarType>::computeKinetic(ForceModel model) const {
        const size_t dof = getDOF();
        Vector<ScalarType> averaged_pos(dof, 0);
        for (size_t i = 0; i < dof; ++i)
            averaged_pos[i] = mean(phasePosX.row(dof + i));

        ScalarType kinetic = repBeta * dof;
        for (size_t replica = 0; replica < getNumReplica(); ++replica) {
            auto col = phasePosX.col(replica);
            auto pos = col.tail(dof);
            MDCellType cell = phaseToCell(replica);
            cell.normalizeCell();
            kinetic += (averaged_pos - pos) * model.force(std::move(cell));
        }
        kinetic /= ScalarType(getNumReplica() * 2);
        return kinetic;
    }

    template<class ScalarType, class PosScalarType>
    ScalarType RPMD<ScalarType, PosScalarType>::estimateTemperature() const {
        return getClassicalKinetic() * (2 / PhyConst<AU>::boltzmannK) / (Dim * getNumParticle() * getNumReplica() * getNumReplica());
    }

    template<class ScalarType, class PosScalarType>
    void RPMD<ScalarType, PosScalarType>::setTemperature(ScalarType temperature) {
        assert(!temperature.isNegative());
        temperatureT = temperature;
        repBeta = temperatureT * PhyConst<AU>::boltzmannK * getNumReplica();
        omegaW = repBeta / PhyConst<AU>::reducedPlanck;
    }

    template<class ScalarType, class PosScalarType>
    void RPMD<ScalarType, PosScalarType>::toNormalRepr(size_t posID) {
        assert(posID < getDOF());
        fft.transform(phasePosX.row(posID));
        auto momentum = buffer.row(0);
        momentum = fft.getKSpace();

        fft.transform(phasePosX.row(posID + getDOF()));
        auto pos = buffer.row(1);
        pos = fft.getKSpace();
    }

    template<class ScalarType, class PosScalarType>
    void RPMD<ScalarType, PosScalarType>::toBeadRepr(size_t posID) {
        assert(posID < getDOF());
        fft.invTransform(buffer.row(0));
        auto momentum = phasePosX.row(posID);
        momentum = fft.getRSpace();

        fft.invTransform(buffer.row(1));
        auto pos = phasePosX.row(posID + getDOF());
        pos = fft.getRSpace();
    }

    template<class ScalarType, class PosScalarType>
    template<class RandomGenerator>
    void RPMD<ScalarType, PosScalarType>::thermostatStep(RandomGenerator& gen, ScalarType deltaT) {
        std::normal_distribution<> dist{};
        const size_t dof = getDOF();
        for (size_t i = 0; i < dof; ++i) {
            const auto mass = cell.getMass(i / Dim);
            const ScalarType factor = sqrt(repBeta * mass * getNumReplica());
            toNormalRepr(i);
            /* Translational mode */ {
                const ScalarType viscosityY = Core::reciprocal(thermostatTime);
                thermostatImpl(0, deltaT, viscosityY, factor, ComplexScalar<ScalarType>(dist(gen)));
            }
            for (size_t j = 1; j < buffer.getColumn(); ++j) {
                const ScalarType phase = M_PI * j / getNumReplica();
                const ScalarType viscosityY = sin(phase) * omegaW;
                const ScalarType normalized_rand = M_SQRT1_2 * dist(gen);
                thermostatImpl(j, deltaT, viscosityY, factor, ComplexScalar<ScalarType>(normalized_rand, normalized_rand));
            }
            toBeadRepr(i);
        }
    }

    template<class ScalarType, class PosScalarType>
    void RPMD<ScalarType, PosScalarType>::thermostatImpl(size_t mode_index, ScalarType deltaT, ScalarType viscosityY, ScalarType factor, ComplexScalar<ScalarType> random) {
        const ScalarType c1 = exp(-viscosityY * deltaT);
        const ScalarType c2 = sqrt(ScalarType(1) - square(c1));
        buffer(0, mode_index) = c1 * buffer(0, mode_index) + factor * c2 * random;
    }

    template<class ScalarType, class PosScalarType>
    void RPMD<ScalarType, PosScalarType>::forceStep(ScalarType deltaT) {
        const size_t dof = getDOF();
        for (size_t replica = 0; replica < getNumReplica(); ++replica) {
            auto phasePos = phasePosX.col(replica);
            auto momentum = phasePos.head(dof);
            momentum += forceBuffer.col(replica).asVector() * deltaT;
        }
    }

    template<class ScalarType, class PosScalarType>
    void RPMD<ScalarType, PosScalarType>::dynamicStep(ScalarType deltaT) {
        using MatrixType = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Element, 2, 2>;
        using VectorType = Vector<ComplexScalar<ScalarType>, 2>;
        const size_t dof = getDOF();

        MatrixType matA(2, 2);
        VectorType temp{};
        for (size_t i = 0; i < dof; ++i) {
            const auto mass = cell.getMass(i / Dim);
            toNormalRepr(i);
            /* Translational mode */ {
                buffer(1, 0) += buffer(0, 0) * deltaT / mass;
            }
            for (size_t j = 1; j < buffer.getColumn(); ++j) {
                auto col = buffer.col(j);
                const ScalarType omegaK = omegaW * sin(ScalarType(M_PI * j / getNumReplica())) * 2;
                const ScalarType factor = ScalarType(mass) * omegaK;
                const ScalarType phase = omegaK * deltaT;
                const ScalarType cosine = cos(phase);
                const ScalarType sine = sin(phase);
                matA(0, 0) = cosine;
                matA(0, 1) = -factor * sine;
                matA(1, 0) = sine / factor;
                matA(1, 1) = cosine;
                temp = matA * col;
                col = temp;
            }
            toBeadRepr(i);
        }
    }

    template<class ScalarType, class PosScalarType>
    void RPMD<ScalarType, PosScalarType>::normalizeCentroid() {
        PositionMatrix centroid = getPos();
        cell.toDirect(centroid);
        size_t index = getDOF();
        for (const auto elem : centroid) {
            const size_t component = index % Dim;
            const size_t atom_start = index - component;
            const int integer = float(elem);
            const Vector<PosScalarType, 3> delta = PosScalarType(integer - elem.isNegative()) * cell.getLattice().row(component).asVector();
            for (size_t i = 0; i < 3; ++i) {
                auto row = phasePosX.row(atom_start + i);
                row -= delta[i];
            }
            ++index;
        }
        assert(checkCentroid());
    }

    template<class ScalarType, class PosScalarType>
    bool RPMD<ScalarType, PosScalarType>::checkCentroid() const {
        constexpr bool success = true;
        PositionMatrix centroid = getPos();
        cell.toDirect(centroid);
        for (auto& elem : centroid)
            if (!(PosScalarType::Zero() <= elem && elem <= PosScalarType::One()))
                return !success;
        return success;
    }
}
