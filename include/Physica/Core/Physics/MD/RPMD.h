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
#include "Physica/Core/Physics/PhyConst.h"
#include "Physica/Core/Math/Statistics/NumCharacter.h"
#include "MDCell.h"
#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] M. Ceriotti, M. Parrinello, T. E. Markland and D. E. Manolopoulos, J. Chem. Phys. 133, 124104 (2010).
     * [2] G. Bussi and M. Parrinello, Phys. Rev. E 75, 056707 (2007).
     */
    template<class ScalarType>
    class RPMD final {
        using PhasePosType = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Vector>;
        using BufferType = DenseMatrix<ComplexScalar<ScalarType>, MatrixOption::Row | MatrixOption::Vector, 2>;
        constexpr static unsigned int Dim = 3;
    private:
        FFT<ScalarType, 1> fft;
        MDCell cell;
        PhasePosType phasePosX;
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
        void step(RandomGenerator& gen, const ForceCalculator& force);
        /* Getters */
        [[nodiscard]] size_t getNumReplica() const noexcept { return phasePosX.getColumn(); }
        [[nodiscard]] size_t getDOF() const noexcept { return Dim * cell.getNumParticle(); }
        [[nodiscard]] typename MDCell::PositionMatrix getPos() const;
        [[nodiscard]] typename MDCell::PositionMatrix getMomentum() const;
        template<class ForceCalculator>
        [[nodiscard]] ScalarType computeKinetic(ForceCalculator force) const;
    private:
        void toNormalRepr(size_t posID);
        void toBeadRepr(size_t posID);
        MDCell phaseToCell(size_t replica) const;
        void cellToPhase(const MDCell& md_cell, size_t replica);
        template<class RandomGenerator>
        void thermostatStep(RandomGenerator& gen);
        template<class ForceCalculator>
        void forceStep(const ForceCalculator& force, unsigned int replica);
        void dynamicStep();
    };

    template<class ScalarType>
    RPMD<ScalarType>::RPMD(MDCell cell_, size_t numReplica, ScalarType temperatureT_, ScalarType thermostatTime_, ScalarType timeStep_)
            : cell(std::move(cell_))
            , temperatureT(std::move(temperatureT_))
            , thermostatTime(std::move(thermostatTime_))
            , timeStep(std::move(timeStep_)) {
        fft = FFT<ScalarType, 1>(numReplica, 1);

        const size_t dof = getDOF();
        phasePosX.resize(2 * dof, numReplica);
        buffer.resize(2, fft.getFreqSize());

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
        repBeta = temperatureT * PhyConst<AU>::boltzmannK * numReplica;
        omegaW = repBeta / PhyConst<AU>::reducedPlanck;
    }

    template<class ScalarType>
    template<class RandomGenerator, class ForceCalculator, class Executor>
    void RPMD<ScalarType>::step(RandomGenerator& gen, const ForceCalculator& force) {
        thermostatStep(gen);
        auto kernel = [&](unsigned int replica) { forceStep(force, replica); };
        Executor::parallel_for(kernel, getNumReplica(), Executor::getNumThread()).wait();
        dynamicStep();
        Executor::parallel_for(kernel, getNumReplica(), Executor::getNumThread()).wait();
        thermostatStep(gen);
    }

    template<class ScalarType>
    typename MDCell::PositionMatrix RPMD<ScalarType>::getPos() const {
        using PositionMatrix = typename MDCell::PositionMatrix;
        using ScalarType_ = typename MDCell::ScalarType;

        PositionMatrix result(cell.getNumParticle(), 3, 0);
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
        using ScalarType_ = typename MDCell::ScalarType;

        PositionMatrix result(cell.getNumParticle(), 3, 0);
        size_t index = 0;
        for (auto& elem : result) {
            elem = ScalarType_(mean(phasePosX.row(index)));
            ++index;
        }
        return result;
    }

    template<class ScalarType>
    template<class ForceCalculator>
    ScalarType RPMD<ScalarType>::computeKinetic(ForceCalculator force) const {
        const size_t dof = getDOF();
        Vector<ScalarType> averaged_pos(dof, 0);
        for (size_t i = 0; i < dof; ++i)
            averaged_pos[i] = mean(phasePosX.row(dof + i));

        ScalarType kinetic = repBeta * dof;
        for (size_t i = 0; i < getNumReplica(); ++i) {
            auto col = phasePosX.col(i);
            auto pos = col.tail(dof);
            kinetic += (averaged_pos - pos) * force(phaseToCell(i));
        }
        kinetic /= ScalarType(2 * getNumReplica());
        return kinetic;
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
    MDCell RPMD<ScalarType>::phaseToCell(size_t replica) const {
        using PositionMatrix = typename MDCell::PositionMatrix;
        using ScalarType_ = typename PositionMatrix::ScalarType;

        PositionMatrix pos(cell.getNumParticle(), 3);
        auto phase = phasePosX.col(replica);
        size_t index = getDOF();
        for (auto& elem : pos) {
            elem = ScalarType_(phase[index]);
            ++index;
        }
        return MDCell(cell.getLattice(), std::move(pos), cell.getMassVec());
    }

    template<class ScalarType>
    void RPMD<ScalarType>::cellToPhase(const MDCell& md_cell, size_t replica) {
        auto phase = phasePosX.col(replica);
        size_t index = getDOF();
        for (auto elem : md_cell.getPos()) {
            phase[index] = elem;
            ++index;
        }
    }

    template<class ScalarType>
    template<class RandomGenerator>
    void RPMD<ScalarType>::thermostatStep(RandomGenerator& gen) {
        std::normal_distribution<> dist{};
        const size_t dof = getDOF();
        for (size_t i = 0; i < dof; ++i) {
            const auto mass = cell.getMass(i / Dim);
            const ScalarType factor = sqrt(repBeta * mass);
            toNormalRepr(i);
            for (size_t j = 0; j < buffer.getColumn(); ++j) {
                ScalarType viscosityY{};
                if (j == 0)
                    viscosityY = reciprocal(thermostatTime);
                else {
                    const ScalarType phase = M_PI * j / getNumReplica();
                    viscosityY = sin(phase) * (omegaW * 2);
                }
                const ScalarType c1 = exp(-viscosityY * (timeStep * 0.5));
                const ScalarType c2 = sqrt(ScalarType(1) - square(c1));
                buffer(0, j) = c1 * buffer(0, j) + factor * c2 * ComplexScalar<ScalarType>(dist(gen), dist(gen));
            }
            toBeadRepr(i);
        }
    }

    template<class ScalarType>
    template<class ForceCalculator>
    void RPMD<ScalarType>::forceStep(const ForceCalculator& force, unsigned int replica) {
        const size_t dof = getDOF();
        auto phasePos = phasePosX.col(replica);
        auto momentum = phasePos.head(dof);
        momentum += force(phaseToCell(replica)) * (timeStep * 0.5);
    }

    template<class ScalarType>
    void RPMD<ScalarType>::dynamicStep() {
        using MatrixType = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Element, 2, 2>;
        using VectorType = Vector<ComplexScalar<ScalarType>, 2>;
        const size_t dof = getDOF();

        MatrixType matA(2, 2);
        VectorType temp{};
        for (size_t i = 0; i < dof; ++i) {
            const auto mass = cell.getMass(i / Dim);
            toNormalRepr(i);
            /* Translational mode */ {
                buffer(1, 0) += buffer(0, 0) * timeStep / mass;
            }
            for (size_t j = 1; j < buffer.getColumn(); ++j) {
                auto col = buffer.col(j);
                const ScalarType omegaK = omegaW * sin(ScalarType(M_PI * j / getNumReplica())) * 2;
                const ScalarType factor = ScalarType(mass) * square(omegaK);
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

        for (size_t replica = 0; replica < getNumReplica(); ++replica) {
            auto cell = phaseToCell(replica);
            cell.checkPeriodic();
            cellToPhase(cell, replica);
        }
    }
}
