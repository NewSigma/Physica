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
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "Physica/Core/Physics/PhyConst.h"

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
        CrystalCell cell;
        PhasePosType phasePosX;
        BufferType buffer;
        ScalarType temperatureT;
        ScalarType thermostatTime;
        ScalarType timeStep;

        ScalarType repBeta;
        ScalarType omegaW;
    public:
        RPMD(CrystalCell cell_, size_t numReplica, ScalarType temperatureT_, ScalarType thermostatTime_, ScalarType timeStep_);
        /* Operations */
        template<class RandomGenerator, class ForceCalculator>
        void step(RandomGenerator gen, ForceCalculator force);
        /* Getters */
        [[nodiscard]] size_t getNumReplica() const noexcept { return phasePosX.getColumn(); }
    private:
        void toNormalRepr(size_t atomId);
        void toBeadRepr(size_t atomId);
        CrystalCell phaseToCell(size_t replica) const;
        template<class RandomGenerator>
        void thermostatStep(RandomGenerator gen);
        template<class ForceCalculator>
        void forceStep(ForceCalculator force);
    };

    template<class ScalarType>
    RPMD<ScalarType>::RPMD(CrystalCell cell_, size_t numReplica, ScalarType temperatureT_, ScalarType thermostatTime_, ScalarType timeStep_)
            : cell(std::move(cell_))
            , temperatureT(std::move(temperatureT_))
            , thermostatTime(std::move(thermostatTime_))
            , timeStep(std::move(timeStep_)) {
        const size_t dof = Dim * cell.getAtomCount();
        fft = FFT<ScalarType, 1>(dof, 1);
        phasePosX.resize(2 * dof, numReplica);
        buffer.resize(2, fft.getFreqSize());

        auto momentum = phasePosX.topRows(dof);
        momentum = ScalarType::Zero();

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

        repBeta = temperatureT * PhyConst<AU>::boltzmannK * numReplica;
        omegaW = repBeta / PhyConst<AU>::reducedPlanck;
    }

    template<class ScalarType>
    template<class RandomGenerator, class ForceCalculator>
    void RPMD<ScalarType>::step(RandomGenerator gen, ForceCalculator force) {
        using MatrixType = DenseMatrix<ComplexScalar<ScalarType>, MatrixOption::Row | MatrixOption::Element, 2, 2>;
        using VectorType = Vector<ComplexScalar<ScalarType>, 2>;
        thermostatStep(std::ref(gen));
        forceStep(std::cref(force));

        MatrixType matA(2, 2);
        VectorType temp{};
        matA(0, 0) = ScalarType(1);
        matA(1, 1) = ScalarType(1);
        for (size_t i = 0; i < cell.getAtomCount(); ++i) {
            const auto atomicNum = cell.getAtomicNumber(i);
            const auto mass = PhyConst<SI>::relativeAtomMass[atomicNum];
            const ScalarType factor = ScalarType(mass) * square(omegaW);

            toNormalRepr(i);
            matA(1, 0) = reciprocal(ScalarType(mass)) * timeStep;
            for (size_t j = 0; j < buffer.getColumn(); ++j) {
                auto col = buffer.col(i);
                const ScalarType phase = 2 * M_PI * j / getNumReplica();
                matA(0, 1) = factor * ComplexScalar<ScalarType>(cos(phase), -sin(phase)) - factor;
                temp = matA * col;
                col = temp;
            }
            toBeadRepr(i);
        }

        forceStep(std::cref(force));
        thermostatStep(std::ref(gen));
    }

    template<class ScalarType>
    void RPMD<ScalarType>::toNormalRepr(size_t atomId) {
        fft.transform(phasePosX.row(atomId));
        auto momentum = buffer.row(0);
        momentum = fft.getFreqs();

        const size_t dof = Dim * cell.getAtomCount();
        fft.transform(phasePosX.row(atomId + dof));
        auto pos = buffer.row(1);
        pos = fft.getFreqs();
    }

    template<class ScalarType>
    void RPMD<ScalarType>::toBeadRepr(size_t atomId) {
        fft.invTransform(buffer.row(0));
        auto momentum = phasePosX.row(atomId);
        momentum = fft.getDatas();

        const size_t dof = Dim * cell.getAtomCount();
        fft.invTransform(buffer.row(1));
        auto pos = phasePosX.row(atomId + dof);
        pos = fft.getDatas();
    }

    template<class ScalarType>
    CrystalCell RPMD<ScalarType>::phaseToCell(size_t replica) const {
        using PositionMatrix = typename CrystalCell::PositionMatrix;
        using ScalarType_ = typename PositionMatrix::ScalarType;

        PositionMatrix pos(cell.getAtomCount(), 3);
        auto phase = phasePosX.col(replica);
        size_t index = Dim * cell.getAtomCount();
        for (auto& elem : pos) {
            elem = ScalarType_(phase[index]);
            ++index;
        }
        return CrystalCell(cell.getLattice(), std::move(pos), cell.getAtomicNumbers());
    }

    template<class ScalarType>
    template<class RandomGenerator>
    void RPMD<ScalarType>::thermostatStep(RandomGenerator gen) {
        std::normal_distribution<> dist{};
        for (size_t i = 0; i < cell.getAtomCount(); ++i) {
            const auto atomicNum = cell.getAtomicNumber(i);
            const auto mass = PhyConst<SI>::relativeAtomMass[atomicNum];
            const ScalarType factor = sqrt(repBeta * mass);
            toNormalRepr(i);
            for (size_t j = 0; j < buffer.getColumn(); ++j) {
                ComplexScalar<ScalarType> viscosityY{};
                if (j == 0)
                    viscosityY = reciprocal(thermostatTime);
                else {
                    const ScalarType phase = 2 * M_PI * j / getNumReplica();
                    viscosityY = sqrt(ComplexScalar<ScalarType>(ScalarType(1) - cos(phase), -sin(phase))) * (omegaW * 2);
                }
                const auto c1 = exp(-viscosityY * (timeStep * 0.5));
                const auto c2 = sqrt(ScalarType(1) - square(c1));
                buffer(0, j) = c1 * buffer(0, j) + factor * c2 * ComplexScalar<ScalarType>(dist(gen.get()), dist(gen.get()));
            }
            toBeadRepr(i);
        }
    }

    template<class ScalarType>
    template<class ForceCalculator>
    void RPMD<ScalarType>::forceStep(ForceCalculator force) {
        const size_t dof = Dim * cell.getAtomCount();
        for (size_t i = 0; i < getNumReplica(); ++i) {
            auto phasePos = phasePosX.col(i);
            auto momentum = phasePos.head(dof);
            momentum -= force(std::move(phaseToCell(i))) * (timeStep * 0.5);
        }
    }
}
