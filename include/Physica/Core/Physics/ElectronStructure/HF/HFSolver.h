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
     */
    template<class BaseSetType>
    class HFSolver {
        using ScalarType = typename Internal::Traits<BaseSetType>::ScalarType;
        using MatrixType = DenseMatrix<ScalarType, DenseMatrixOption::Column | DenseMatrixOption::Vector>;
    private:
        const Molecular<ScalarType>& molecular;
        MatrixType singleHamilton;
        MatrixType overlap;
        Utils::Array<BaseSetType> baseSet;
        ScalarType selfConsistentEnergy;
        size_t electronCount;
    public:
        HFSolver(const Molecular<ScalarType>& m, size_t electronCount_, size_t baseSetSize);
        HFSolver(const HFSolver&) = delete;
        HFSolver(HFSolver&&) noexcept = delete;
        ~HFSolver() = default;
        /* Operators */
        HFSolver& operator=(const HFSolver& base) = delete;
        HFSolver& operator=(HFSolver&& base) noexcept = delete;
        /* Operations */
        bool compute(const ScalarType& criteria, size_t maxIte);
        /* Getters */
        [[nodiscard]] Utils::Array<BaseSetType>& getBaseSet() noexcept { return baseSet; }
        [[nodiscard]] const Utils::Array<BaseSetType>& getBaseSet() const noexcept { return baseSet; }
        [[nodiscard]] size_t getBaseSetSize() const noexcept { return baseSet.getLength(); }
        [[nodiscard]] ScalarType getSelfConsistentEnergy() const noexcept { return selfConsistentEnergy; }
        [[nodiscard]] ScalarType getTotalEnergy() const noexcept { return selfConsistentEnergy + molecular.getNuclearRepulsionEnergy(); }
    private:
        void formSingleHamilton();
        void formOverlapMatrix();
        void formDensityMatrix(MatrixType& density, const MatrixType& wave_func);
        void formCoulombMatrix(MatrixType& __restrict fock, const MatrixType& __restrict density);
        template<class VectorType>
        void updateSelfConsistentEnergy(const VectorType& eigenvalues, const Utils::Array<size_t>& sortedEigenvalues, const MatrixType& waves);
        void sortEigenvalues(const typename EigenSolver<MatrixType>::EigenvalueVector& eigenvalues, Utils::Array<size_t>& indexToSort) const;
    };

    template<class BaseSetType>
    HFSolver<BaseSetType>::HFSolver(const Molecular<ScalarType>& m, size_t electronCount_, size_t baseSetSize)
            : molecular(m)
            , singleHamilton(baseSetSize, baseSetSize)
            , overlap(baseSetSize, baseSetSize)
            , baseSet(baseSetSize)
            , selfConsistentEnergy()
            , electronCount(electronCount_) {
        assert(electronCount <= baseSetSize);
    }
    /**
     * Perform self-consistant computation
     * 
     * \return true if converged, false otherwise
     */
    template<class BaseSetType>
    bool HFSolver<BaseSetType>::compute(const ScalarType& criteria, size_t maxIte) {
        assert(criteria.isPositive());

        const size_t baseSetSize = getBaseSetSize();
        formSingleHamilton();
        formOverlapMatrix();

        const MatrixType cholesky = Cholesky(overlap);
        const MatrixType inv_cholesky = cholesky.inverse();

        MatrixType densityMat = MatrixType::Zeros(baseSetSize);
        MatrixType fock = singleHamilton;
        MatrixType waves = MatrixType(baseSetSize, electronCount);
        Utils::Array<size_t> sortedEigenvalues = Utils::Array<size_t>(getBaseSetSize());

        size_t iteration = 0;
        do {
            const MatrixType modifiedFock = (inv_cholesky * fock).compute() * inv_cholesky.transpose();
            const EigenSolver<MatrixType> solver(modifiedFock, true);
            
            const auto& eigenvalues = solver.getEigenvalues();
            sortEigenvalues(eigenvalues, sortedEigenvalues);
            auto eigenvectors = solver.getEigenvectors();
            for (size_t i = 0; i < electronCount; ++i) {
                auto wave = waves.col(i);
                wave.asVector() = (inv_cholesky.transpose() * toRealVector(eigenvectors.col(sortedEigenvalues[i])).moveToColMatrix()).compute().col(0);
            }
            formDensityMatrix(densityMat, waves);
            formCoulombMatrix(fock, densityMat);
            // Get ground state energy
            const ScalarType oldSelfConsistentEnergy = selfConsistentEnergy;
            updateSelfConsistentEnergy(eigenvalues, sortedEigenvalues, waves);
            const ScalarType delta = abs(oldSelfConsistentEnergy - selfConsistentEnergy);
            std::cout << iteration << ' ' << oldSelfConsistentEnergy << std::endl;
            // Check convergence
            if (delta < criteria)
                return true;
            if ((++iteration) == maxIte)
                return false;
            // Prepare for next iteration
            fock += singleHamilton;
        } while(true);
        return true;
    }

    template<class BaseSetType>
    void HFSolver<BaseSetType>::formSingleHamilton() {
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
    void HFSolver<BaseSetType>::formOverlapMatrix() {
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
    void HFSolver<BaseSetType>::formDensityMatrix(MatrixType& density, const MatrixType& wave_func) {
        const size_t baseSetSize = getBaseSetSize();
        for (size_t i = 0; i < baseSetSize; ++i) {
            size_t j = 0;
            for (; j < i; ++j)
                density(j, i) = density(i, j);

            for (; j < baseSetSize; ++j) {
                ScalarType temp = ScalarType::Zero();
                for (size_t k = 0; k < electronCount; ++k) {
                    auto wave = wave_func.col(k);
                    temp += wave[i] * wave[j];
                }
                density(j, i) = ScalarType::Two() * temp;
            }
        }
    }

    template<class BaseSetType>
    void HFSolver<BaseSetType>::formCoulombMatrix(MatrixType& __restrict fock, const MatrixType& __restrict density) {
        const size_t size = getBaseSetSize();
        for (size_t p = 0; p < size; ++p) {
            for (size_t q = 0; q < size; ++q) {
                ScalarType temp = ScalarType::Zero();
                for (size_t r = 0; r < size; ++r) {
                    for (size_t s = 0; s < size; ++s) {
                        const ScalarType coulomb = BaseSetType::electronRepulsion(baseSet[p], baseSet[r], baseSet[q], baseSet[s]);
                        const ScalarType exchange = BaseSetType::electronRepulsion(baseSet[p], baseSet[r], baseSet[s], baseSet[q]);
                        temp += density(s, r) * (coulomb - ScalarType(0.5) * exchange);
                    }
                }
                fock(q, p) = temp;
            }
        }
    }

    template<class BaseSetType>
    template<class VectorType>
    void HFSolver<BaseSetType>::updateSelfConsistentEnergy(const VectorType& eigenvalues,
                                                           const Utils::Array<size_t>& sortedEigenvalues,
                                                           const MatrixType& waves) {
        selfConsistentEnergy = ScalarType::Zero();
        for (size_t i = 0; i < electronCount; ++i) {
            selfConsistentEnergy += eigenvalues[sortedEigenvalues[i]].getReal();
            auto wave = waves.col(i);
            selfConsistentEnergy += ((waves.transpose() * singleHamilton).compute() * wave).calc(0, 0);
        }
        selfConsistentEnergy *= ScalarType(0.5);
    }
    /**
     * Get the first \param electronCount lowest eigenvalues and save their indexes to array \param index,
     * the eigenvalues are in ascending order.
     * 
     * \param index
     * A array whose length is \param electronCount
     */
    template<class BaseSetType>
    void HFSolver<BaseSetType>::sortEigenvalues(const typename EigenSolver<MatrixType>::EigenvalueVector& eigenvalues,
                                                Utils::Array<size_t>& indexToSort) const {
        auto arrayToSort = toRealVector(eigenvalues);
        for (size_t i = 0; i < getBaseSetSize(); ++i)
            indexToSort[i] = i;

        for (size_t i = 0; i < electronCount; ++i) {
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
