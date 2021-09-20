/*
 * Copyright 2021 WeiBo He.
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

#include "Physica/Core/Physics/Molecular.h"
#include "Physica/Core/Physics/ElectronStructure/ElectronConfig.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/Cholesky.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/InverseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Transpose.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/EigenSolver.h"

namespace Physica::Core::Physics {
    namespace Internal {
        template<class T> class Traits;
    }
    /**
     * Reference:
     * [1] Jos Thijssen. Computational Physics[M].London: Cambridge university press, 2013:43-88
     * [2] Larsen A, Poulsen R S. Applied Hartree-Fock methods.
     */
    template<class BaseSetType>
    class RHFSolver {
        using ScalarType = typename Internal::Traits<BaseSetType>::ScalarType;
        using MatrixType = DenseMatrix<ScalarType, DenseMatrixOption::Column | DenseMatrixOption::Vector>;
        /**
         * Infer next wave function from the \param DIISMultiplicity wave functions before, refer to [2]
         */
        constexpr static size_t DIISMultiplicity = 3;
        static_assert(DIISMultiplicity >= 3, "DIISMultiplicity less than three makes no sence");
        using WaveGroup = Utils::Array<MatrixType, DIISMultiplicity>;
        using DeltaWaveGroup = Utils::Array<MatrixType, DIISMultiplicity - 1>;
        using DIISMatrix = DenseMatrix<ScalarType, DenseMatrixOption::Column | DenseMatrixOption::Element, DIISMultiplicity, DIISMultiplicity>;
    private:
        const Molecular<ScalarType>& molecular;
        ElectronConfig electronConfig;
        size_t numOccupiedOrbit;
        MatrixType singleHamilton;
        MatrixType overlap;
        Utils::Array<BaseSetType> baseSet;
        ScalarType selfConsistentEnergy;
        EigenSolver<MatrixType> eigenSolver;
        size_t iteration;
    public:
        RHFSolver(const Molecular<ScalarType>& m, const ElectronConfig& electronConfig_, size_t baseSetSize);
        RHFSolver(const RHFSolver&) = delete;
        RHFSolver(RHFSolver&&) noexcept = delete;
        ~RHFSolver() = default;
        /* Operators */
        RHFSolver& operator=(const RHFSolver& base) = delete;
        RHFSolver& operator=(RHFSolver&& base) noexcept = delete;
        /* Operations */
        bool compute(const ScalarType& criteria, size_t maxIte);
        /* Getters */
        [[nodiscard]] Utils::Array<BaseSetType>& getBaseSet() noexcept { return baseSet; }
        [[nodiscard]] const Utils::Array<BaseSetType>& getBaseSet() const noexcept { return baseSet; }
        [[nodiscard]] size_t getBaseSetSize() const noexcept { return baseSet.getLength(); }
        [[nodiscard]] ScalarType getSelfConsistentEnergy() const noexcept { return selfConsistentEnergy; }
        [[nodiscard]] ScalarType getTotalEnergy() const noexcept { return selfConsistentEnergy + molecular.getNuclearRepulsionEnergy(); }
        [[nodiscard]] size_t getIteration() const noexcept { return iteration; }
    private:
        void formSingleHamilton();
        void formOverlapMatrix();
        void formDensityMatrix(MatrixType& __restrict electronDensity,
                               MatrixType& __restrict sameSpinElectronDensity,
                               const MatrixType& __restrict wave_func);
        void formCoulombMatrix(MatrixType& __restrict fock,
                               const MatrixType& __restrict electronDensity,
                               const MatrixType& __restrict sameSpinElectronDensity);
        void updateWaves(const MatrixType& inv_cholesky,
                         WaveGroup& waveGroup,
                         DeltaWaveGroup& deltaWaveGroup,
                         const Utils::Array<size_t>& sortedEigenvalues);
        [[nodiscard]] static ScalarType waveMatrixDot(const MatrixType& m1, const MatrixType& m2);
        void updateSelfConsistentEnergy(const Utils::Array<size_t>& sortedEigenvalues, const MatrixType& waveGroup);
        void sortEigenvalues(Utils::Array<size_t>& indexToSort) const;
    };

    template<class BaseSetType>
    RHFSolver<BaseSetType>::RHFSolver(const Molecular<ScalarType>& m, const ElectronConfig& electronConfig_, size_t baseSetSize)
            : molecular(m)
            , electronConfig(electronConfig_)
            , numOccupiedOrbit(0)
            , singleHamilton(baseSetSize, baseSetSize)
            , overlap(baseSetSize, baseSetSize)
            , baseSet(baseSetSize)
            , selfConsistentEnergy()
            , eigenSolver(baseSetSize)
            , iteration(0) {
        numOccupiedOrbit = electronConfig.getNumOccupiedOrbit();
        assert(numOccupiedOrbit <= baseSetSize);
    }
    /**
     * Perform self-consistant computation
     * 
     * \return true if converged, false otherwise
     */
    template<class BaseSetType>
    bool RHFSolver<BaseSetType>::compute(const ScalarType& criteria, size_t maxIte) {
        assert(criteria.isPositive());

        const size_t baseSetSize = getBaseSetSize();
        formSingleHamilton();
        formOverlapMatrix();

        const MatrixType cholesky = Cholesky(overlap);
        const MatrixType inv_cholesky = cholesky.inverse();

        MatrixType electronDensity = MatrixType::Zeros(baseSetSize);
        MatrixType sameSpinElectronDensity = MatrixType::Zeros(baseSetSize);
        MatrixType fock = singleHamilton;
        WaveGroup waveGroup = WaveGroup(DIISMultiplicity, MatrixType::Zeros(baseSetSize, numOccupiedOrbit));
        DeltaWaveGroup deltaWaveGroup = DeltaWaveGroup(DIISMultiplicity - 1, MatrixType::Zeros(baseSetSize, numOccupiedOrbit));
        Utils::Array<size_t> sortedEigenvalues = Utils::Array<size_t>(getBaseSetSize());

        iteration = 0;
        do {
            const MatrixType modifiedFock = (inv_cholesky * fock).compute() * inv_cholesky.transpose();
            eigenSolver.compute(modifiedFock, true);

            sortEigenvalues(sortedEigenvalues);
            updateWaves(inv_cholesky, waveGroup, deltaWaveGroup, sortedEigenvalues);
            // Get ground state energy
            const ScalarType oldSelfConsistentEnergy = selfConsistentEnergy;
            updateSelfConsistentEnergy(sortedEigenvalues, *waveGroup.rbegin());
            const ScalarType delta = abs(oldSelfConsistentEnergy - selfConsistentEnergy);
            // Check convergence
            if (delta < criteria)
                return true;
            if ((++iteration) == maxIte)
                return false;
            // Prepare for next iteration
            formDensityMatrix(electronDensity, sameSpinElectronDensity, *waveGroup.rbegin());
            formCoulombMatrix(fock, electronDensity, sameSpinElectronDensity);
            fock += singleHamilton;
        } while(true);
        return true;
    }

    template<class BaseSetType>
    void RHFSolver<BaseSetType>::formSingleHamilton() {
        const size_t baseSetSize = getBaseSetSize();
        for (size_t i = 0; i < baseSetSize; ++i) {
            size_t j = 0;
            for (; j < i; ++j)
                singleHamilton(j, i) = singleHamilton(i, j);

            for (; j < baseSetSize; ++j) {
                ScalarType temp = ScalarType::Zero();
                for (size_t k = 0; k < molecular.getAtomCount(); ++k)
                    temp -= BaseSetType::nuclearAttraction(baseSet[i], baseSet[j], molecular.getAtom(k).getVector())
                            * ScalarType(molecular.getAtomicNumber(k));
                singleHamilton(j, i) = BaseSetType::kinetic(baseSet[i], baseSet[j]) + temp;
            }
        }
    }

    template<class BaseSetType>
    void RHFSolver<BaseSetType>::formOverlapMatrix() {
        const size_t baseSetSize = getBaseSetSize();
        for (size_t i = 0; i < baseSetSize; ++i) {
            size_t j = 0;
            for (; j < i; ++j)
                overlap(j, i) = overlap(i, j);

            for (; j < baseSetSize; ++j)
                overlap(j, i) = BaseSetType::overlap(baseSet[i], baseSet[j]);
        }
    }

    template<class BaseSetType>
    void RHFSolver<BaseSetType>::formDensityMatrix(MatrixType& __restrict electronDensity,
                                                   MatrixType& __restrict sameSpinElectronDensity,
                                                   const MatrixType& __restrict wave_func) {
        const size_t baseSetSize = getBaseSetSize();
        for (size_t i = 0; i < baseSetSize; ++i) {
            size_t j = 0;
            for (; j < i; ++j) {
                electronDensity(j, i) = electronDensity(i, j);
                sameSpinElectronDensity(j, i) = sameSpinElectronDensity(i, j);
            }

            for (; j < baseSetSize; ++j) {
                ScalarType temp1 = ScalarType::Zero();
                ScalarType temp2 = ScalarType::Zero();
                for (size_t k = 0; k < numOccupiedOrbit; ++k) {
                    const size_t orbitPos = electronConfig.getOccupiedOrbitPos(k);
                    const auto orbitState = electronConfig.getOrbitState(orbitPos);
                    assert(orbitState != ElectronConfig::NoOccupacy);
                    const bool isSingleOccupacy = orbitState == ElectronConfig::SingleOccupacy;
                    auto wave = wave_func.col(k);
                    const ScalarType dot = wave[i] * wave[j];
                    temp1 += isSingleOccupacy ? dot : ScalarType::Two() * dot;
                    temp2 += dot;
                }
                electronDensity(j, i) = temp1;
                sameSpinElectronDensity(j, i) = temp2;
            }
        }
    }

    template<class BaseSetType>
    void RHFSolver<BaseSetType>::formCoulombMatrix(MatrixType& __restrict fock,
                                                   const MatrixType& __restrict electronDensity,
                                                   const MatrixType& __restrict sameSpinElectronDensity) {
        const size_t size = getBaseSetSize();
        for (size_t p = 0; p < size; ++p) {
            for (size_t q = 0; q < size; ++q) {
                ScalarType temp = ScalarType::Zero();
                for (size_t r = 0; r < size; ++r) {
                    for (size_t s = 0; s < size; ++s) {
                        const ScalarType coulomb = BaseSetType::electronRepulsion(baseSet[p], baseSet[r], baseSet[q], baseSet[s]);
                        const ScalarType exchange = BaseSetType::electronRepulsion(baseSet[p], baseSet[r], baseSet[s], baseSet[q]);
                        temp += electronDensity(s, r) * coulomb - sameSpinElectronDensity(s, r) * exchange;
                    }
                }
                fock(q, p) = temp;
            }
        }
    }

    template<class BaseSetType>
    void RHFSolver<BaseSetType>::updateWaves(const MatrixType& inv_cholesky,
                                             WaveGroup& waveGroup,
                                             DeltaWaveGroup& deltaWaveGroup,
                                             const Utils::Array<size_t>& sortedEigenvalues) {
        for (size_t i = 0; i < waveGroup.getLength() - 2; ++i) {
            swap(waveGroup[i], waveGroup[i + 1]);
            swap(deltaWaveGroup[i], deltaWaveGroup[i + 1]);
        }
        swap(waveGroup[waveGroup.getLength() - 2], waveGroup[waveGroup.getLength() - 1]);

        auto& newWave = *waveGroup.rbegin();
        const auto& eigenvectors = eigenSolver.getRawEigenvectors();
        for (size_t i = 0; i < numOccupiedOrbit; ++i) {
            auto eigenState = newWave.col(i);
            const size_t orbitPos = electronConfig.getOccupiedOrbitPos(i);
            const size_t solutionPos = sortedEigenvalues[orbitPos];
            eigenState.asVector() = (inv_cholesky.transpose() * eigenvectors.col(solutionPos)).compute().col(0);
        }
        deltaWaveGroup[waveGroup.getLength() - 2] = newWave - waveGroup[waveGroup.getLength() - 2];

        const bool readyForDIIS = iteration >= DIISMultiplicity - 1 && iteration % 6 == 0; // Constant 6 comes from Reference: Improved SCF Convergence Acceleration[J]. Journal of Computational Chemistry, 1982, 3(4):556-560.
        if (readyForDIIS) {
            Vector<ScalarType, DIISMultiplicity> x{};
            /* Solve linear equation */ {
                DIISMatrix A = DIISMatrix(DIISMultiplicity, DIISMultiplicity, -ScalarType::One());
                A(DIISMultiplicity - 1, DIISMultiplicity - 1) = ScalarType::Zero();
                for (size_t i = 0; i < A.getColumn() - 1; ++i)
                    for (size_t j = 0; j < A.getRow() - 1; ++j)
                        A(i, j) = waveMatrixDot(deltaWaveGroup[i], deltaWaveGroup[j]);
                Vector<ScalarType, DIISMultiplicity> b = Vector<ScalarType, DIISMultiplicity>(DIISMultiplicity, ScalarType::Zero());
                b[DIISMultiplicity - 1] = -ScalarType::One();
                const DIISMatrix inv_A = A.inverse();
                x = (inv_A * b.moveToColMatrix()).compute().col(0);
            }
            newWave = ScalarType::Zero();
            for (size_t i = 0; i < DIISMultiplicity - 1; ++i)
                newWave += waveGroup[i] * x[i];

            for (size_t i = 0; i < numOccupiedOrbit; ++i) {
                auto col = newWave.col(i);
                auto norm = ((col.transpose() * overlap).compute() * col).calc(0, 0);
                col.asVector() *= reciprocal(sqrt(norm));
            }
        }
    }

    template<class BaseSetType>
    typename RHFSolver<BaseSetType>::ScalarType RHFSolver<BaseSetType>::waveMatrixDot(const MatrixType& m1, const MatrixType& m2) {
        ScalarType result = ScalarType::Zero();
        assert(m1.getRow() == m2.getRow());
        assert(m1.getColumn() == m2.getColumn());
        for (size_t i = 0; i < m1.getColumn(); ++i) {
            auto col1 = m1.col(i);
            auto col2 = m2.col(i);
            result += col1.asVector() * col2.asVector();
        }
        return result;
    }

    template<class BaseSetType>
    void RHFSolver<BaseSetType>::updateSelfConsistentEnergy(const Utils::Array<size_t>& sortedEigenvalues,
                                                            const MatrixType& waveGroup) {
        const auto& eigenvalues = eigenSolver.getEigenvalues();
        selfConsistentEnergy = ScalarType::Zero();
        for (size_t i = 0; i < waveGroup.getColumn(); ++i) {
            const size_t orbitPos = electronConfig.getOccupiedOrbitPos(i);
            ScalarType temp = eigenvalues[sortedEigenvalues[orbitPos]].getReal();
            auto wave = waveGroup.col(i);
            temp += ((wave.transpose() * singleHamilton).compute() * wave).calc(0, 0);
            const auto orbitState = electronConfig.getOrbitState(orbitPos);
            assert(orbitState != ElectronConfig::NoOccupacy);
            const bool isSingleOccupacy = orbitState == ElectronConfig::SingleOccupacy;
            selfConsistentEnergy += isSingleOccupacy ? temp : (ScalarType::Two() * temp);
        }
        selfConsistentEnergy *= ScalarType(0.5);
    }
    /**
     * Get the first \param orbitCount lowest eigenvalues and save their indexes to array \param index,
     * the eigenvalues are in ascending order.
     * 
     * \param index
     * A array whose length is \param orbitCount
     */
    template<class BaseSetType>
    void RHFSolver<BaseSetType>::sortEigenvalues(Utils::Array<size_t>& indexToSort) const {
        const auto& eigenvalues = eigenSolver.getEigenvalues();
        auto arrayToSort = toRealVector(eigenvalues);
        for (size_t i = 0; i < getBaseSetSize(); ++i)
            indexToSort[i] = i;

        for (size_t i = 0; i < getBaseSetSize(); ++i) {
            size_t indexOfToInsert = i;
            for (size_t j = i + 1; j < getBaseSetSize(); ++j) {
                if (arrayToSort[indexOfToInsert] > arrayToSort[j])
                    indexOfToInsert = j;
            }
            std::swap(arrayToSort[i], arrayToSort[indexOfToInsert]);
            std::swap(indexToSort[i], indexToSort[indexOfToInsert]);
            assert(eigenvalues[indexToSort[i]].getReal() >= eigenvalues[indexToSort[i == 0 ? 0 : i - 1]].getReal());
        }
    }
}
