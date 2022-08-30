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
     */
    template<class ScalarType>
    class RPMD final : private MDCell {
        using PhasePosType = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Vector>;
        using ForceMatrix = DenseMatrix<ScalarType, MatrixOption::Column | MatrixOption::Vector>;
        using BufferType = DenseMatrix<ComplexScalar<ScalarType>, MatrixOption::Row | MatrixOption::Vector, 2>;
        constexpr static unsigned int Dim = 3;
    private:
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
        RPMD(MDCell cell_, size_t numReplica, ScalarType temperatureT_, ScalarType thermostatTime_, ScalarType timeStep_);
        /* Operations */
        template<class RandomGenerator, class ForceCalculator, class Executor = Parallel::SequentialExecutor>
        void nvt_step(RandomGenerator& gen, const ForceCalculator& force);
        template<class ForceCalculator, class Executor = Parallel::SequentialExecutor>
        void nve_step(const ForceCalculator& force);
        template<class RandomGenerator, class ForceCalculator, class Executor = Parallel::SequentialExecutor>
        void nvt_step_for(ScalarType duration, RandomGenerator& gen, const ForceCalculator& force);
        template<class ForceCalculator, class Executor = Parallel::SequentialExecutor>
        void nve_step_for(ScalarType duration, const ForceCalculator& force);
        template<class RandomGenerator>
        void initMomentum(RandomGenerator& gen);
        void scaleVelocity();
        /* Getters */
        [[nodiscard]] const PhasePosType& getPhasePos() const noexcept { return phasePosX; }
        [[nodiscard]] PhasePosType& getPhasePos() noexcept { return phasePosX; }
        [[nodiscard]] size_t getNumReplica() const noexcept { return phasePosX.getColumn(); }
        [[nodiscard]] MDCell phaseToCell(size_t replica) const;
        [[nodiscard]] size_t getDOF() const noexcept { return Dim * MDCell::getNumParticle(); }
        [[nodiscard]] ScalarType getTemperature() const noexcept { return temperatureT; }
        [[nodiscard]] typename MDCell::PositionMatrix getPos() const;
        [[nodiscard]] typename MDCell::PositionMatrix getMomentum() const;
        [[nodiscard]] const ForceMatrix& getForce() const noexcept { return forceBuffer; }

        [[nodiscard]] ScalarType getClassicalKinetic() const;
        template<class ForceModel>
        [[nodiscard]] ScalarType getClassicalPotentialEnergy(ForceModel model) const;
        template<class ForceModel>
        [[nodiscard]] ScalarType getClassicalInternalEnergy(ForceModel model) const;
        template<class ForceCalculator>
        [[nodiscard]] ScalarType computeKinetic(ForceCalculator force) const;
        [[nodiscard]] ScalarType estimateTemperature() const;
        /* Setters */
        void setTemperature(ScalarType temperature);
    private:
        void toNormalRepr(size_t posID);
        void toBeadRepr(size_t posID);
        template<class RandomGenerator>
        void thermostatStep(RandomGenerator& gen);
        void thermostatImpl(size_t mode_index, ScalarType viscosityY, ScalarType factor, ComplexScalar<ScalarType> random);
        void forceStep();
        void dynamicStep();
        void normalizeCentroid();
        bool checkCentroid() const;
    };

    template<class ScalarType>
    RPMD<ScalarType>::RPMD(MDCell cell_, size_t numReplica, ScalarType temperatureT_, ScalarType thermostatTime_, ScalarType timeStep_)
            : MDCell(std::move(cell_))
            , thermostatTime(std::move(thermostatTime_))
            , timeStep(std::move(timeStep_)) {
        fft = FFT<ScalarType, 1>(numReplica, 1);

        const size_t dof = getDOF();
        phasePosX.resize(2 * dof, numReplica);
        forceBuffer.resize(dof, numReplica);
        buffer.resize(2, fft.getFreqSize());

        auto momentum = phasePosX.topRows(dof);
        momentum = ScalarType::Zero();

        /* Fill pos */ {
            size_t index = dof;
            for (auto elem : MDCell::getPos()) {
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

    template<class ScalarType>
    template<class RandomGenerator, class ForceCalculator, class Executor>
    void RPMD<ScalarType>::nvt_step(RandomGenerator& gen, const ForceCalculator& force) {
        thermostatStep(gen);
        nve_step<ForceCalculator, Executor>(force);
        thermostatStep(gen);
    }

    template<class ScalarType>
    template<class ForceCalculator, class Executor>
    void RPMD<ScalarType>::nve_step(const ForceCalculator& force) {
        auto kernel = [&](unsigned int replica) {
            MDCell cell = phaseToCell(replica);
            cell.normalizeCell();
            auto saveTo = forceBuffer.col(replica);
            saveTo = force(std::move(cell));
        };
        Executor::parallel_for(kernel, getNumReplica(), Executor::getNumThread()).wait();
        forceStep();
        dynamicStep();
        normalizeCentroid();
        Executor::parallel_for(kernel, getNumReplica(), Executor::getNumThread()).wait();
        forceStep();
    }

    template<class ScalarType>
    template<class RandomGenerator, class ForceCalculator, class Executor>
    void RPMD<ScalarType>::nvt_step_for(ScalarType duration, RandomGenerator& gen, const ForceCalculator& force) {
        uint64_t step = double(duration / timeStep) + 0.5;
        for (uint64_t _ = 0; _ < step; ++_)
            nvt_step<RandomGenerator, ForceCalculator, Executor>(gen, force);
    }

    template<class ScalarType>
    template<class ForceCalculator, class Executor>
    void RPMD<ScalarType>::nve_step_for(ScalarType duration, const ForceCalculator& force) {
        uint64_t step = double(duration / timeStep) + 0.5;
        for (uint64_t _ = 0; _ < step; ++_)
            nve_step<ForceCalculator, Executor>(force);
    }

    template<class ScalarType>
    template<class RandomGenerator>
    void RPMD<ScalarType>::initMomentum(RandomGenerator& gen) {
        std::normal_distribution<> dist{};
        const size_t dof = getDOF();
        Vector<ScalarType, 3> driftMomentum(3, 0);
        for (size_t i = 0; i < dof; ++i) {
            const auto mass = MDCell::getMass(i / Dim);
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

    template<class ScalarType>
    void RPMD<ScalarType>::scaleVelocity() {
        const ScalarType temperatureNow = estimateTemperature();
        assert(temperatureNow.isPositive());
        const size_t dof = getDOF();
        const ScalarType factor = sqrt(temperatureT / temperatureNow);
        auto momentum = phasePosX.topRows(dof);
        momentum *= factor;
    }

    template<class ScalarType>
    MDCell RPMD<ScalarType>::phaseToCell(size_t replica) const {
        using PositionMatrix = typename MDCell::PositionMatrix;
        using ScalarType_ = typename PositionMatrix::ScalarType;

        PositionMatrix pos(MDCell::getNumParticle(), 3);
        auto phase = phasePosX.col(replica);
        size_t index = getDOF();
        for (auto& elem : pos) {
            elem = ScalarType_(phase[index]);
            ++index;
        }
        return MDCell(MDCell::getLattice(), std::move(pos), MDCell::getMassVec());
    }

    template<class ScalarType>
    typename MDCell::PositionMatrix RPMD<ScalarType>::getPos() const {
        using PositionMatrix = typename MDCell::PositionMatrix;
        using ScalarType_ = typename PositionMatrix::ScalarType;

        PositionMatrix result(MDCell::getNumParticle(), 3);
        const size_t dof = getDOF();
        size_t index = dof;
        for (auto& elem : result) {
            elem = ScalarType_(mean(phasePosX.row(index)));
            ++index;
        }
        return result;
    }

    template<class ScalarType>
    typename MDCell::PositionMatrix RPMD<ScalarType>::getMomentum() const {
        using PositionMatrix = typename MDCell::PositionMatrix;
        using ScalarType_ = typename PositionMatrix::ScalarType;

        PositionMatrix result(MDCell::getNumParticle(), 3, 0);
        size_t index = 0;
        for (auto& elem : result) {
            elem = ScalarType_(mean(phasePosX.row(index)));
            ++index;
        }
        return result;
    }

    template<class ScalarType>
    ScalarType RPMD<ScalarType>::getClassicalKinetic() const {
        const size_t dof = getDOF();
        ScalarType classical_kinetic = 0;
        for (size_t i = 0; i < dof; ++i) {
            const auto mass = MDCell::getMass(i / Dim);
            auto p = phasePosX.row(i);
            classical_kinetic += square(p.asVector()).sum() / (mass * 2);
        }
        return classical_kinetic;
    }

    template<class ScalarType>
    template<class ForceModel>
    ScalarType RPMD<ScalarType>::getClassicalPotentialEnergy(ForceModel model) const {
        ScalarType result = 0;
        for (size_t i = 0; i < getNumReplica(); ++i)
            result += model.potentialEnergy(phaseToCell(i));
        return result;
    }

    template<class ScalarType>
    template<class ForceModel>
    ScalarType RPMD<ScalarType>::getClassicalInternalEnergy(ForceModel model) const {
        return getClassicalKinetic() + getClassicalPotentialEnergy<ForceModel>(model);
    }

    template<class ScalarType>
    template<class ForceCalculator>
    ScalarType RPMD<ScalarType>::computeKinetic(ForceCalculator force) const {
        const size_t dof = getDOF();
        Vector<ScalarType> averaged_pos(dof, 0);
        for (size_t i = 0; i < dof; ++i)
            averaged_pos[i] = mean(phasePosX.row(dof + i));

        ScalarType kinetic = repBeta * dof;
        for (size_t replica = 0; replica < getNumReplica(); ++replica) {
            auto col = phasePosX.col(replica);
            auto pos = col.tail(dof);
            MDCell cell = phaseToCell(replica);
            cell.normalizeCell();
            kinetic += (averaged_pos - pos) * force(std::move(cell));
        }
        kinetic /= ScalarType(getNumReplica() * 2);
        return kinetic;
    }

    template<class ScalarType>
    ScalarType RPMD<ScalarType>::estimateTemperature() const {
        return getClassicalKinetic() * (2 / PhyConst<AU>::boltzmannK) / (Dim * getNumParticle() * getNumReplica() * getNumReplica());
    }

    template<class ScalarType>
    void RPMD<ScalarType>::setTemperature(ScalarType temperature) {
        assert(!temperature.isNegative());
        temperatureT = temperature;
        repBeta = temperatureT * PhyConst<AU>::boltzmannK * getNumReplica();
        omegaW = repBeta / PhyConst<AU>::reducedPlanck;
    }

    template<class ScalarType>
    void RPMD<ScalarType>::toNormalRepr(size_t posID) {
        assert(posID < getDOF());
        fft.transform(phasePosX.row(posID));
        auto momentum = buffer.row(0);
        momentum = fft.getFreqs();

        fft.transform(phasePosX.row(posID + getDOF()));
        auto pos = buffer.row(1);
        pos = fft.getFreqs();
    }

    template<class ScalarType>
    void RPMD<ScalarType>::toBeadRepr(size_t posID) {
        assert(posID < getDOF());
        fft.invTransform(buffer.row(0));
        auto momentum = phasePosX.row(posID);
        momentum = fft.getDatas();

        fft.invTransform(buffer.row(1));
        auto pos = phasePosX.row(posID + getDOF());
        pos = fft.getDatas();
    }

    template<class ScalarType>
    template<class RandomGenerator>
    void RPMD<ScalarType>::thermostatStep(RandomGenerator& gen) {
        std::normal_distribution<> dist{};
        const size_t dof = getDOF();
        for (size_t i = 0; i < dof; ++i) {
            const auto mass = MDCell::getMass(i / Dim);
            const ScalarType factor = sqrt(repBeta * mass * getNumReplica());
            toNormalRepr(i);
            /* Translational mode */ {
                const ScalarType viscosityY = Core::reciprocal(thermostatTime);
                thermostatImpl(0, viscosityY, factor, ComplexScalar<ScalarType>(dist(gen)));
            }
            for (size_t j = 1; j < buffer.getColumn(); ++j) {
                const ScalarType phase = M_PI * j / getNumReplica();
                const ScalarType viscosityY = sin(phase) * omegaW;
                const ScalarType normalized_rand = M_SQRT1_2 * dist(gen);
                thermostatImpl(j, viscosityY, factor, ComplexScalar<ScalarType>(normalized_rand, normalized_rand));
            }
            toBeadRepr(i);
        }
    }

    template<class ScalarType>
    void RPMD<ScalarType>::thermostatImpl(size_t mode_index, ScalarType viscosityY, ScalarType factor, ComplexScalar<ScalarType> random) {
        const ScalarType c1 = exp(-viscosityY * (timeStep * 0.5));
        const ScalarType c2 = sqrt(ScalarType(1) - square(c1));
        buffer(0, mode_index) = c1 * buffer(0, mode_index) + factor * c2 * random;
    }

    template<class ScalarType>
    void RPMD<ScalarType>::forceStep() {
        const size_t dof = getDOF();
        for (size_t replica = 0; replica < getNumReplica(); ++replica) {
            auto phasePos = phasePosX.col(replica);
            auto momentum = phasePos.head(dof);
            momentum += forceBuffer.col(replica).asVector() * (timeStep * 0.5);
        }
    }

    template<class ScalarType>
    void RPMD<ScalarType>::dynamicStep() {
        using MatrixType = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Element, 2, 2>;
        using VectorType = Vector<ComplexScalar<ScalarType>, 2>;
        const size_t dof = getDOF();

        MatrixType matA(2, 2);
        VectorType temp{};
        for (size_t i = 0; i < dof; ++i) {
            const auto mass = MDCell::getMass(i / Dim);
            toNormalRepr(i);
            /* Translational mode */ {
                buffer(1, 0) += buffer(0, 0) * timeStep / mass;
            }
            for (size_t j = 1; j < buffer.getColumn(); ++j) {
                auto col = buffer.col(j);
                const ScalarType omegaK = omegaW * sin(ScalarType(M_PI * j / getNumReplica())) * 2;
                const ScalarType factor = ScalarType(mass) * omegaK;
                const ScalarType phase = omegaK * timeStep;
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

    template<class ScalarType>
    void RPMD<ScalarType>::normalizeCentroid() {
        using ScalarType_ = typename PositionMatrix::ScalarType;
        PositionMatrix centroid = getPos();
        MDCell::toDirect(centroid);
        size_t index = getDOF();
        for (const auto elem : centroid) {
            const size_t component = index % Dim;
            const size_t atom_start = index - component;
            const int integer = float(elem);
            const Vector<ScalarType_, 3> delta = ScalarType_(integer - elem.isNegative()) * MDCell::getLattice().row(component).asVector();
            for (size_t i = 0; i < 3; ++i) {
                auto row = phasePosX.row(atom_start + i);
                row -= delta[i];
            }
            ++index;
        }
        assert(checkCentroid());
    }

    template<class ScalarType>
    bool RPMD<ScalarType>::checkCentroid() const {
        using ScalarType_ = typename PositionMatrix::ScalarType;
        constexpr bool success = true;
        PositionMatrix centroid = getPos();
        MDCell::toDirect(centroid);
        for (auto& elem : centroid)
            if (!(ScalarType_::Zero() <= elem && elem <= ScalarType_::One()))
                return !success;
        return success;
    }
}
