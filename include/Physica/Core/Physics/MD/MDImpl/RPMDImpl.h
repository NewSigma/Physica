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

#include "Physica/Core/Math/Calculus/ODE/SRK2.h"

namespace Physica::Core {
    template<class ScalarType, class PosScalarType, unsigned int Dim>
    RPMD<ScalarType, PosScalarType, Dim>::RPMD(MDCellType cell_,
                                               size_t numReplica,
                                               size_t numContract,
                                               ScalarType temperatureT_,
                                               ScalarType thermostatTime_,
                                               ScalarType timeStep_)
            : cell(std::move(cell_))
            , fft(numReplica, 1)
            , fftContract(numContract, 1)
            , thermostatTime(std::move(thermostatTime_))
            , timeStep(std::move(timeStep_)) {
        assert(0 < numContract && numContract <= numReplica);

        const size_t dof = getDOF();
        phaseMatrix.resize(2 * dof, numReplica);
        forceBuffer.resize(dof, numReplica);
        if (isContractEnabled()) {
            posContract.resize(dof, numContract);
            forceContract.resize(dof, numContract);
        }
        buffer.resize(2, fft.getKSpaceSize());

        auto momentum = phaseMatrix.topRows(dof);
        momentum = ScalarType::Zero();

        /* Fill pos */ {
            size_t index = dof;
            for (auto elem : cell.getPos()) {
                phaseMatrix(index, 0) = elem;
                ++index;
            }
            for (size_t i = 1; i < getNumReplica(); ++i) {
                auto col = phaseMatrix.col(i);
                auto pos = col.tail(dof);
                pos = phaseMatrix.col(0).tail(dof);
            }
        }
        setTemperature(temperatureT_);
        checkParam();
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class ForceModel, class Executor>
    void RPMD<ScalarType, PosScalarType, Dim>::updateForce(const ForceModel& model) {
        if (isContractEnabled()) {
            auto kernel_short = [&](unsigned int replica) {
                MDCellType cell = phaseToCell(replica);
                cell.normalize();
                auto saveTo = forceBuffer.col(replica);
                saveTo = model.template force_short<Executor>(std::move(cell));
            };
            auto future_short = Executor::parallel_for(kernel_short, getNumReplica(), Executor::getNumThread());

            contract();
            auto kernel_long = [&](unsigned int contract) {
                MDCellType cell = contractToCell(contract);
                cell.normalize();
                auto saveTo = forceContract.col(contract);
                saveTo = model.template force_long<Executor>(std::move(cell));
            };
            auto future_long = Executor::parallel_for(kernel_long, getNumContract(), Executor::getNumThread());

            Executor::auto_wait(future_long);
            decontract();
            Executor::auto_wait(future_short);
        }
        else {
            auto kernel = [&](unsigned int replica) {
                MDCellType cell = phaseToCell(replica);
                cell.normalize();
                auto saveTo = forceBuffer.col(replica);
                saveTo = model.template force<Executor>(std::move(cell));
            };
            auto future = Executor::parallel_for(kernel, getNumReplica(), Executor::getNumThread());
            Executor::auto_wait(future);
        }
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class RandomGenerator, class ForceModel, class Executor>
    void RPMD<ScalarType, PosScalarType, Dim>::nvt_step(RandomGenerator& gen, const ForceModel& model) {
        forceStep(timeStep * 0.5);
        dynamicStep(timeStep * 0.5);
        thermostatStep(gen, timeStep);
        dynamicStep(timeStep * 0.5);
        updateForce<ForceModel, Executor>(model);
        forceStep(timeStep * 0.5);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class ForceModel, class Executor>
    void RPMD<ScalarType, PosScalarType, Dim>::nve_step(const ForceModel& model) {
        forceStep(timeStep * 0.5);
        dynamicStep(timeStep);
        updateForce<ForceModel, Executor>(model);
        forceStep(timeStep * 0.5);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class RandomGenerator, class Barostat, class ForceModel, class Executor>
    void RPMD<ScalarType, PosScalarType, Dim>::npt_step(RandomGenerator& gen, Barostat& barostat, const ForceModel& model) {
        barostat.forceStep(*this, timeStep * 0.5);
        npt_dynamicStep(barostat, timeStep * 0.5);
        thermostatStep(gen, timeStep);
        npt_dynamicStep(barostat, timeStep * 0.5);
        updateForce<ForceModel, Executor>(model);
        barostat.forceStep(*this, timeStep * 0.5);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class RandomGenerator, class ForceModel, class Executor>
    void RPMD<ScalarType, PosScalarType, Dim>::nvt_step_for(ScalarType duration, RandomGenerator& gen, const ForceModel& model) {
        uint64_t step = double(duration / timeStep) + 0.5;
        for (uint64_t _ = 0; _ < step; ++_)
            nvt_step<RandomGenerator, ForceModel, Executor>(gen, model);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class ForceModel, class Executor>
    void RPMD<ScalarType, PosScalarType, Dim>::nve_step_for(ScalarType duration, const ForceModel& force) {
        uint64_t step = double(duration / timeStep) + 0.5;
        for (uint64_t _ = 0; _ < step; ++_)
            nve_step<ForceModel, Executor>(force);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class RandomGenerator, class Barostat, class ForceModel, class Executor>
    void RPMD<ScalarType, PosScalarType, Dim>::npt_step_for(
            ScalarType duration,
            RandomGenerator& gen,
            Barostat& barostat,
            const ForceModel& model) {
        uint64_t step = double(duration / timeStep) + 0.5;
        for (uint64_t _ = 0; _ < step; ++_)
            npt_step<RandomGenerator, Barostat, ForceModel, Executor>(gen, barostat, model);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class RandomGenerator>
    void RPMD<ScalarType, PosScalarType, Dim>::initMomentum(RandomGenerator& gen) {
        std::normal_distribution<> dist{};
        const size_t dof = getDOF();
        Vector<ScalarType, Dim> driftMomentum(Dim, 0);
        for (size_t i = 0; i < dof; ++i) {
            const auto mass = cell.getMass(i / Dim);
            const size_t direction = i % Dim;
            const ScalarType factor = sqrt(repBeta * mass);
            for (size_t j = 0; j < getNumReplica(); ++j) {
                const ScalarType temp = factor * dist(gen);
                phaseMatrix(i, j) = temp;
                driftMomentum[direction] += temp;
            }
        }
        driftMomentum *= Core::reciprocal(ScalarType(getNumParticle() * getNumReplica()));

        for (size_t i = 0; i < dof; ++i) {
            auto row = phaseMatrix.row(i);
            row -= driftMomentum[i % Dim];
        }
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::removeDrift() {
        const size_t dof = getDOF();
        Vector<ScalarType, Dim> driftMomentum(Dim, 0);
        for (size_t i = 0; i < dof; ++i) {
            const size_t direction = i % Dim;
            for (size_t j = 0; j < getNumReplica(); ++j)
                driftMomentum[direction] += phaseMatrix(i, j);
        }
        driftMomentum *= Core::reciprocal(ScalarType(getNumParticle() * getNumReplica()));

        for (size_t i = 0; i < dof; ++i) {
            auto row = phaseMatrix.row(i);
            row -= driftMomentum[i % Dim];
        }
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::scaleVelocity() {
        const ScalarType temperatureNow = calcTemperature();
        assert(temperatureNow.isPositive());
        const size_t dof = getDOF();
        const ScalarType factor = sqrt(temperatureT / temperatureNow);
        auto momentum = phaseMatrix.topRows(dof);
        momentum *= factor;
    }
    /**
     * Carrying out this function every several steps may stable the simulation.
     */
    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::normalizeCentroid() {
        PositionMatrix centroid = makeCentroidPos();
        cell.toDirect(centroid);
        size_t index = getDOF();
        for (const auto elem : centroid) {
            const size_t component = index % Dim;
            const size_t atom_start = index - component;
            const int integer = float(elem);
            const Vector<PosScalarType, Dim> delta = PosScalarType(integer - elem.isNegative()) * cell.getLattice().row(component).asVector();
            for (size_t i = 0; i < Dim; ++i) {
                auto row = phaseMatrix.row(atom_start + i);
                row -= delta[i];
            }
            ++index;
        }
        assert(checkCentroid());
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    typename RPMD<ScalarType, PosScalarType, Dim>::MDCellType RPMD<ScalarType, PosScalarType, Dim>::phaseToCell(size_t replica) const {
        PositionMatrix pos(getNumParticle(), Dim);
        auto phase = phaseMatrix.col(replica);
        size_t index = getDOF();
        for (auto& elem : pos) {
            elem = PosScalarType(phase[index]);
            ++index;
        }
        return MDCellType(cell.getLattice(), std::move(pos), cell.getMassVec());
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    typename RPMD<ScalarType, PosScalarType, Dim>::MDCellType RPMD<ScalarType, PosScalarType, Dim>::contractToCell(size_t contract) const {
        PositionMatrix pos(getNumParticle(), Dim);
        auto phase = posContract.col(contract);
        size_t index = 0;
        for (auto& elem : pos) {
            elem = PosScalarType(phase[index]);
            ++index;
        }
        return MDCellType(cell.getLattice(), std::move(pos), cell.getMassVec());
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::checkParam() const {
        const ScalarType cycle = ScalarType(2 * M_PI) / omegaW;
        bool isSmallEnough = timeStep < cycle / ScalarType(4);
        if (!isSmallEnough)
            throw std::invalid_argument("[Error]: Time step is too large");
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    typename RPMD<ScalarType, PosScalarType, Dim>::PositionMatrix RPMD<ScalarType, PosScalarType, Dim>::makeCentroidPos() const {
        PositionMatrix result(getNumParticle(), Dim);
        const size_t dof = getDOF();
        size_t index = dof;
        for (auto& elem : result) {
            elem = PosScalarType(mean(phaseMatrix.row(index)));
            ++index;
        }
        return result;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    typename RPMD<ScalarType, PosScalarType, Dim>::MDCellType RPMD<ScalarType, PosScalarType, Dim>::makeAverageCell() const {
        return MDCellType(getLattice(), makeCentroidPos(), cell.getMassVec());
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    typename RPMD<ScalarType, PosScalarType, Dim>::PositionMatrix RPMD<ScalarType, PosScalarType, Dim>::makeCentroidMomentum() const {
        PositionMatrix result(getNumParticle(), Dim, 0);
        size_t index = 0;
        for (auto& elem : result) {
            elem = PosScalarType(mean(phaseMatrix.row(index)));
            ++index;
        }
        return result;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    typename RPMD<ScalarType, PosScalarType, Dim>::PositionMatrix RPMD<ScalarType, PosScalarType, Dim>::getMomentum(size_t replica) const {
        PositionMatrix result(getNumParticle(), Dim, 0);
        size_t index = 0;
        for (auto& elem : result) {
            elem = PosScalarType(phaseMatrix(index, replica));
            ++index;
        }
        return result;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    ScalarType RPMD<ScalarType, PosScalarType, Dim>::getClassicalKinetic() const {
        const size_t dof = getDOF();
        ScalarType classical_kinetic = 0;
        for (size_t i = 0; i < dof; ++i) {
            const auto mass = cell.getMass(i / Dim);
            auto p = phaseMatrix.row(i);
            classical_kinetic += square(p.asVector()).sum() / (mass * 2);
        }
        return classical_kinetic;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class ForceModel>
    ScalarType RPMD<ScalarType, PosScalarType, Dim>::getClassicalPotentialEnergy(const ForceModel& model) const {
        ScalarType result = 0;
        for (size_t i = 0; i < getNumReplica(); ++i)
            result += model.potentialEnergy(phaseToCell(i));
        return result;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    ScalarType RPMD<ScalarType, PosScalarType, Dim>::getClassicalElastic() const {
        const size_t dof = getDOF();
        auto pos = phaseMatrix.bottomRows(dof);
        ScalarType result = 0;
        for (size_t i = 0; i < dof; ++i) {
            const ScalarType mass = cell.getMass(i / Dim);
            for (size_t j = 0; j < getNumReplica(); ++j)
                result += mass * square(omegaW * (pos(i, j) - pos(i, (j + 1) % getNumReplica()))) * 0.5;
        }
        return result;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class ForceModel>
    ScalarType RPMD<ScalarType, PosScalarType, Dim>::getClassicalInternalEnergy(const ForceModel& model) const {
        return getClassicalKinetic() + getClassicalPotentialEnergy<ForceModel>(model) + getClassicalElastic();
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    ScalarType RPMD<ScalarType, PosScalarType, Dim>::calcKinetic() const {
        const size_t dof = getDOF();
        Vector<ScalarType> averaged_pos(dof, 0);
        for (size_t i = 0; i < dof; ++i)
            averaged_pos[i] = mean(phaseMatrix.row(dof + i));

        ScalarType kinetic = repBeta * dof;
        for (size_t replica = 0; replica < getNumReplica(); ++replica) {
            auto phase = phaseMatrix.col(replica);
            auto pos = phase.tail(dof);
            kinetic += (averaged_pos - pos) * forceBuffer.col(replica);
        }
        kinetic /= ScalarType(getNumReplica() * 2);
        return kinetic;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class ForceModel, class Executor>
    ScalarType RPMD<ScalarType, PosScalarType, Dim>::calcPotential(const ForceModel& model) const {
        Vector<ScalarType> temp(getNumReplica());
        auto kernel = [this, model, &temp](unsigned int replica) {
            temp[replica] = model.potentialEnergy(phaseToCell(replica));
        };
        Executor::parallel_for(kernel, getNumReplica(), Executor::getNumThread()).wait();
        return mean(temp);
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    ScalarType RPMD<ScalarType, PosScalarType, Dim>::calcTemperature() const {
        return square(makeCentroidMomentum()).sum() * (1 / (Dim * PhyConst<AU>::boltzmannK)) / cell.getMassVec().sum();
    }
    /**
     * The function has size effect, extend the cell shall ease the problem, refer to [1].
     * 
     * Reference:
     * [1] M. J. Louwerse and E. J. Baerends, Chem. Phys. Lett. 421, 138 (2006).
     */
    template<class ScalarType, class PosScalarType, unsigned int Dim>
    typename RPMD<ScalarType, PosScalarType, Dim>::LatticeMatrix RPMD<ScalarType, PosScalarType, Dim>::makeStress() const {
        assert(getNumReplica() == 1);
        const size_t dof = getDOF();
        LatticeMatrix stress(Dim, Dim, 0);
        for (size_t replica = 0; replica < 1; ++replica) {
            const auto col = phaseMatrix.col(replica);
            const auto momentum = col.head(dof);
            const auto pos = col.tail(dof);
            const auto force = forceBuffer.col(replica);
            for (size_t i = 0; i < getNumParticle(); ++i) {
                const size_t from = i * getDim();
                const size_t to = from + getDim();
                const ScalarType factor = reciprocal(getMassVec()[i]);
                const auto momentum1 = momentum.segment(from, to);
                const auto pos1 = pos.segment(from, to);
                const auto force1 = force.segment(from, to);
                stress += (factor * momentum1) * momentum1.transpose() - pos1 * force1.transpose();
            }
        }
        stress *= reciprocal(getVolume());
        return stress;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::setTemperature(ScalarType temperature) {
        assert(!temperature.isNegative());
        temperatureT = temperature;
        repBeta = temperatureT * PhyConst<AU>::boltzmannK * getNumReplica();
        omegaW = repBeta / PhyConst<AU>::reducedPlanck;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::toNormalRepr(size_t posID) {
        assert(posID < getDOF());
        fft.transform(phaseMatrix.row(posID));
        auto momentum = buffer.row(0);
        momentum = fft.getKSpace();

        fft.transform(phaseMatrix.row(posID + getDOF()));
        auto pos = buffer.row(1);
        pos = fft.getKSpace();
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::toBeadRepr(size_t posID) {
        assert(posID < getDOF());
        fft.invTransform(buffer.row(0));
        auto momentum = phaseMatrix.row(posID);
        momentum = fft.getRSpace();

        fft.invTransform(buffer.row(1));
        auto pos = phaseMatrix.row(posID + getDOF());
        pos = fft.getRSpace();
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::toContractBeadRepr(size_t posID) {
        assert(posID < getDOF());
        auto row = buffer.row(1);
        auto head = row.head(fftContract.getKSpaceSize());
        fftContract.invTransform(head);
        auto pos = posContract.row(posID);
        pos = fftContract.getRSpace() * (ScalarType(getNumContract()) / ScalarType(getNumReplica()));
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::forceToNormRepr(size_t posID) {
        assert(posID < getDOF());
        fftContract.transform(forceContract.row(posID));
        auto row = buffer.row(0);
        row = ScalarType(0);
        auto head = row.head(fftContract.getKSpaceSize());
        head = fftContract.getKSpace();
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::forceToBeadRepr(size_t posID) {
        assert(posID < getDOF());
        fft.invTransform(buffer.row(0));
        auto f = forceBuffer.row(posID);
        f += fft.getRSpace() * (ScalarType(getNumReplica()) / ScalarType(getNumContract()));
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::contract() {
        const size_t dof = getDOF();
        for (size_t i = 0; i < dof; ++i) {
            toNormalRepr(i);
            toContractBeadRepr(i);
        }
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::decontract() {
        const size_t dof = getDOF();
        for (size_t i = 0; i < dof; ++i) {
            forceToNormRepr(i);
            forceToBeadRepr(i);
        }
    }
    /**
     * Reference:
     * [1] G, Bussi, D. Donadio and M. Parrinello, J. Chem. Phys. 126, 014101 (2007).
     */
    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class RandomGenerator>
    void RPMD<ScalarType, PosScalarType, Dim>::thermostatStep(RandomGenerator& gen, ScalarType deltaT) {
        std::normal_distribution<> dist{};
        const size_t dof = getDOF();
        ScalarType factor_translational;
        {
            using VectorType = Vector<ScalarType, 1>;
            [[maybe_unused]] ScalarType _ = 0;
            const ScalarType nowT = calcTemperature();
            VectorType sol{nowT};
            SRK2<ScalarType, 1>::step(timeStep, _, sol,
                                        [this]([[maybe_unused]] ScalarType x, VectorType sol) -> VectorType {
                                            return {(temperatureT - sol[0]) / thermostatTime};
                                        },
                                        [this, &gen]([[maybe_unused]] ScalarType x, VectorType sol) -> VectorType {
                                            std::normal_distribution dist{};
                                            return {sqrt((temperatureT * sol[0]) / (thermostatTime * getDOF())) * 2 * dist(gen)};
                                        });
            factor_translational = sqrt(temperatureT / sol[0]);
        }
        for (size_t i = 0; i < dof; ++i) {
            const auto mass = cell.getMass(i / Dim);
            toNormalRepr(i);
            buffer(0, 0) *= factor_translational;

            const ScalarType factor = sqrt(repBeta * mass * getNumReplica());
            for (size_t j = 1; j < buffer.getColumn(); ++j) {
                const ScalarType phase = M_PI * j / getNumReplica();
                const ScalarType viscosityY = sin(phase) * omegaW;
                const ScalarType normalized_rand = M_SQRT1_2 * dist(gen);
                thermostatImpl(j, deltaT, viscosityY, factor, ComplexScalar<ScalarType>(normalized_rand, normalized_rand));
            }
            toBeadRepr(i);
        }
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::thermostatImpl(size_t mode_index, ScalarType deltaT, ScalarType viscosityY, ScalarType factor, ComplexScalar<ScalarType> random) {
        const ScalarType c1 = exp(-viscosityY * deltaT);
        const ScalarType c2 = sqrt(ScalarType(1) - square(c1));
        buffer(0, mode_index) = c1 * buffer(0, mode_index) + factor * c2 * random;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::forceStep(ScalarType deltaT) {
        const size_t dof = getDOF();
        for (size_t replica = 0; replica < getNumReplica(); ++replica) {
            auto col = phaseMatrix.col(replica);
            auto momentum = col.head(dof);
            momentum += forceBuffer.col(replica).asVector() * deltaT;
        }
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void RPMD<ScalarType, PosScalarType, Dim>::dynamicStep(ScalarType deltaT) {
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

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class Barostat>
    void RPMD<ScalarType, PosScalarType, Dim>::npt_dynamicStep(Barostat& barostat, ScalarType deltaT) {
        dynamicStep(deltaT);
        LatticeMatrix lattice = getLattice() + barostat.getLatticeMomentum() * (deltaT / barostat.getLatticeMass());
        cell.setLattice(std::move(lattice));
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    bool RPMD<ScalarType, PosScalarType, Dim>::checkCentroid() const {
        constexpr bool success = true;
        PositionMatrix centroid = makeCentroidPos();
        cell.toDirect(centroid);
        for (auto& elem : centroid)
            if (!(PosScalarType::Zero() <= elem && elem <= PosScalarType::One()))
                return !success;
        return success;
    }
}
