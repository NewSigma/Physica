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
        void formFockMatrix(MatrixType& __restrict fock, const MatrixType& __restrict density);
        template<class VectorType>
        void updateSelfConsistentEnergy(const VectorType& eigenvalues, const size_t* lowEigenvalueIndex, const MatrixType& electronDencity);
        void getLowestEigenvalues(const typename EigenSolver<MatrixType>::EigenvalueVector& eigenvalues, size_t* index) const;
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

        MatrixType cholesky = Cholesky(overlap);
        MatrixType inv_cholesky = cholesky.inverse();

        MatrixType densityMat = MatrixType::Zeros(baseSetSize);
        MatrixType fock = MatrixType::Zeros(baseSetSize, baseSetSize);
        MatrixType waves = MatrixType(baseSetSize, electronCount);
        size_t lowEigenvalueIndex[electronCount];
        Vector<ScalarType> temp{};

        size_t iteration = 0;
        do {
            fock += singleHamilton;
            const MatrixType modifiedFock = (inv_cholesky * fock).compute() * inv_cholesky.transpose();
            const EigenSolver<MatrixType> solver(modifiedFock, true);
            // Get ground state energy
            const auto& eigenvalues = solver.getEigenvalues();
            getLowestEigenvalues(eigenvalues, lowEigenvalueIndex);

            const ScalarType oldSelfConsistentEnergy = selfConsistentEnergy;
            updateSelfConsistentEnergy(eigenvalues, lowEigenvalueIndex, densityMat);
            const ScalarType delta = abs(oldSelfConsistentEnergy - selfConsistentEnergy);
            // Check convergence
            if (delta < criteria)
                return true;
            if ((++iteration) == maxIte)
                return false;
            // Prepare for next iteration
            auto eigenvectors = solver.getEigenvectors();
            for (size_t i = 0; i < electronCount; ++i) {
                temp = (inv_cholesky * toRealVector(eigenvectors.col(lowEigenvalueIndex[i])).moveToColMatrix()).compute().col(0);
                nomalizeWave(temp);
                auto col = waves.col(i);
                temp.assignTo(col);
            }
            formDensityMatrix(densityMat, waves);
            formFockMatrix(fock, densityMat);
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
        const size_t maxMajor = density.getMaxMajor();
        const size_t maxMinor = density.getMaxMinor();
        for (size_t i = 0; i < maxMajor; ++i) {
            auto col_i = wave_func.col(i);
            size_t j = 0;
            for (; j < i; ++j)
                density.getElementFromMajorMinor(i, j) = density.getElementFromMajorMinor(j, i);

            for (; j < maxMinor; ++j) {
                auto col_j = wave_func.col(j);
                ScalarType temp = ScalarType::Zero();
                for (size_t k = 0; k < electronCount; ++k)
                    temp += col_i[k] * col_j[k];
                density.getElementFromMajorMinor(i, j) = ScalarType::Two() * temp;
            }
        }
    }

    template<class BaseSetType>
    void HFSolver<BaseSetType>::formFockMatrix(MatrixType& __restrict fock, const MatrixType& __restrict density) {
        const size_t size = getBaseSetSize();
        for (size_t p = 0; p < size; ++p) {
            for (size_t q = 0; q < size; ++q) {
                ScalarType temp = ScalarType::Zero();
                for (size_t r; r < size; ++r) {
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
    void HFSolver<BaseSetType>::updateSelfConsistentEnergy(const VectorType& eigenvalues, const size_t* lowEigenvalueIndex, const MatrixType& electronDencity) {
        const size_t baseSetSize = getBaseSetSize();
        selfConsistentEnergy = ScalarType::Zero();
        for (size_t i = 0; i < electronCount; ++i)
            selfConsistentEnergy += eigenvalues[lowEigenvalueIndex[i]].getReal();

        selfConsistentEnergy *= ScalarType(0.5);
        for (size_t i = 0; i < baseSetSize; ++i)
            for (size_t j = i; j < baseSetSize; ++j)
                selfConsistentEnergy += singleHamilton(j, i) * electronDencity(j, i);
    }

    template<class BaseSetType>
    void HFSolver<BaseSetType>::getLowestEigenvalues(const typename EigenSolver<MatrixType>::EigenvalueVector& eigenvalues, size_t* index) const {
        size_t i = 0;
        size_t max_index = 0;
        for (; i < electronCount; ++i) {
            assert(eigenvalues[i].getImag().isZero());
            index[i] = i;
            max_index = eigenvalues[i].getReal() > eigenvalues[max_index].getReal() ? i : max_index;
        }

        for (; i < getBaseSetSize(); ++i) {
            assert(eigenvalues[i].getImag().isZero());
            const auto& value = eigenvalues[i].getReal();
            if (value < eigenvalues[i].getReal()) {
                index[max_index] = i;
                for (size_t j = 0; j < electronCount; ++j)
                    max_index = eigenvalues[j].getReal() > eigenvalues[max_index].getReal() ? j : max_index;
            }
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
