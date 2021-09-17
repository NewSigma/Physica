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
        void updateSelfConsistentEnergy(const VectorType& eigenvalues, const Utils::Array<size_t>& sortedEigenvalues, const MatrixType& electronDencity);
        void sortEigenvalues(const typename EigenSolver<MatrixType>::EigenvalueVector& eigenvalues, Utils::Array<size_t>& indexToSort) const;
        void nomalizeWave(Vector<ScalarType>& eigenvector);
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
        MatrixType fock = MatrixType::Zeros(baseSetSize, baseSetSize);
        MatrixType older_waves = MatrixType(baseSetSize, electronCount);
        MatrixType old_waves = MatrixType(baseSetSize, electronCount);
        MatrixType waves = MatrixType(baseSetSize, electronCount);
        Utils::Array<size_t> sortedEigenvalues = Utils::Array<size_t>(electronCount);
        Vector<ScalarType> temp{};
        Vector<ScalarType> temp1{};

        size_t iteration = 0;
        do {
            fock += singleHamilton;
            const MatrixType modifiedFock = (inv_cholesky * fock).compute() * inv_cholesky.transpose();
            const EigenSolver<MatrixType> solver(modifiedFock, true);
            // Get ground state energy
            const auto& eigenvalues = solver.getEigenvalues();
            sortEigenvalues(eigenvalues, sortedEigenvalues);

            const ScalarType oldSelfConsistentEnergy = selfConsistentEnergy;
            updateSelfConsistentEnergy(eigenvalues, sortedEigenvalues, densityMat);
            const ScalarType delta = abs(oldSelfConsistentEnergy - selfConsistentEnergy);
            //std::cout << iteration << ' ' << oldSelfConsistentEnergy << std::endl;
            // Check convergence
            if (delta < criteria)
                return true;
            if ((++iteration) == maxIte)
                return false;
            // Prepare for next iteration
            auto eigenvectors = solver.getEigenvectors();
            for (size_t i = 0; i < electronCount; ++i) {
                temp = (inv_cholesky.transpose() * toRealVector(eigenvectors.col(sortedEigenvalues[i])).moveToColMatrix()).compute().col(0);
                nomalizeWave(temp);
                auto older_wave = older_waves.col(i);
                auto old_wave = old_waves.col(i);
                auto wave = waves.col(i);
                if (iteration > 3)
                    wave.asVector() = divide(multiply(temp, old_wave) - square(older_wave), temp + older_wave - older_wave.asVector() * ScalarType::Two());
                else
                    wave = temp;
            }
            formDensityMatrix(densityMat, waves);
            older_waves.swap(old_waves);
            old_waves.swap(waves);
            formCoulombMatrix(fock, densityMat);
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
                                                           const MatrixType& electronDencity) {
        const size_t baseSetSize = getBaseSetSize();
        selfConsistentEnergy = ScalarType::Zero();
        for (size_t i = 0; i < electronCount; ++i)
            selfConsistentEnergy += eigenvalues[sortedEigenvalues[i]].getReal();

        for (size_t i = 0; i < baseSetSize; ++i)
            for (size_t j = i; j < baseSetSize; ++j)
                selfConsistentEnergy += singleHamilton(j, i) * electronDencity(j, i);
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
        indexToSort[0] = 0;
        for (size_t i = 1; i < electronCount; ++i) {
            const ScalarType& toInsert = eigenvalues[i].getReal();
            size_t insertTo = 0;
            for (; insertTo < i; ++insertTo) {
                if (toInsert < eigenvalues[insertTo].getReal())
                    break;
            }
            size_t temp = i;
            for (size_t j = insertTo; j <= i; ++j)
                std::swap(temp, indexToSort[j]);
        }
    }

    template<class BaseSetType>
    void HFSolver<BaseSetType>::nomalizeWave(Vector<ScalarType>& eigenvector) {
        auto mat = eigenvector.moveToColMatrix();
        ScalarType norm = ((mat.transpose() * overlap).compute() * mat).calc(0, 0);
        mat *= reciprocal(norm);
        eigenvector = mat.col(0);
    }
}
