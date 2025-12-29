/*
 * Copyright 2021-2025 Weibo He.
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

#include "Physica/Core/Physics/Molecular/Molecular.h"
#include "Physica/Core/Physics/ElectronicStructure/ElectronConfig.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/Cholesky.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/EigenSolver.h"
#include "Physica/Core/Math/Optimization/QuadraticProgramming/QuadraticProgramming.h"

namespace Physica {
    /**
     * Reference:
     * [1] J. H. Thijssen. Computational Physics[M]. London: Cambridge University Press, 2013:43-88
     * [2] Larsen A, Poulsen R S. Applied Hartree-Fock methods.
     * [3] J. Chem. Phys. 116(19), 8255-8261 (2002); https://doi.org/10.1063/1.1470195
     */
    template<class BaseSetType>
    class RHFSolver {
        using This = RHFSolver<BaseSetType>;
        using T = Traits<BaseSetType>::ScalarType;

        constexpr static size_t EDIISBufferSize = 3; //Refer EDIIS from [3]
        constexpr static size_t DIISBufferSize = 3; //Refer DIIS from [2]
        constexpr static size_t MatrixBufferSize = std::max(EDIISBufferSize, DIISBufferSize);
        static_assert(DIISBufferSize >= 3, "DIISBufferSize less than three makes no sence");
        using EDIISBuffer = Array<MatrixND<T>, EDIISBufferSize>;
        using DIISBuffer = Array<MatrixND<T>, DIISBufferSize - 1>;
        using MatrixBuffer = Array<MatrixND<T>, MatrixBufferSize>;
        using DIISMatrix = DenseMatrix<T, MatrixOption::Col, DIISBufferSize, DIISBufferSize>;
    private:
        const Molecular<T>& molecular;
        ElectronConfig electronConfig;
        size_t numOccupiedOrbit;
        DenseSymmMatrix<T, Dynamic> singleHamilton;
        MatrixND<T> overlap;
        Array<BaseSetType> baseSet;
        T selfConsistentEnergy;
        MatrixND<T> wave;
        EigenSolver<T> eigenSolver;
        size_t iteration = 0;
    public:
        RHFSolver(const Molecular<T>& m, const ElectronConfig& electronConfig_, size_t baseSetSize);
        RHFSolver(const RHFSolver&) = delete;
        RHFSolver(RHFSolver&&) noexcept = delete;
        ~RHFSolver() = default;
        /* Operators */
        RHFSolver& operator=(const RHFSolver& base) = delete;
        RHFSolver& operator=(RHFSolver&& base) noexcept = delete;
        /* Operations */
        bool compute(const T& criteria, size_t maxIte);
        /* Getters */
        [[nodiscard]] auto&& getBaseSet(this auto&&) noexcept;
        [[nodiscard]] size_t getBaseSetSize() const noexcept { return baseSet.getLength(); }
        [[nodiscard]] T getSelfConsistentEnergy() const noexcept { return selfConsistentEnergy; }
        [[nodiscard]] T getTotalEnergy() const noexcept { return selfConsistentEnergy + molecular.getNuclearRepulsionEnergy(); }
        [[nodiscard]] const auto& getWave() const noexcept { return wave; }
        [[nodiscard]] size_t getIteration() const noexcept { return iteration; }
        /* Setters */
        void setInitialWave(const MatrixND<T>& initialWave) { wave = initialWave; }
    private:
        void formSingleHamilton();
        void formOverlapMatrix();
        void formDensityMatrix(EDIISBuffer& densityMatrices,
                               MatrixND<T>& sameSpinElectronDensity);
        void formFockMatrix(MatrixBuffer& fockMatrices,
                            const MatrixND<T>& electronDensity,
                            const MatrixND<T>& sameSpinElectronDensity);
        void preDIIS(const MatrixBuffer& fockMatrices,
                     DIISBuffer& errorMatrices,
                     const MatrixND<T>& electronDensity,
                     const MatrixND<T>& inv_cholesky,
                     DIISMatrix& DIISMat);
        void EDIISInterpolation(MatrixBuffer& fockMatrices,
                                EDIISBuffer& densityMatrices,
                                const DenseVector<T, EDIISBufferSize>& energyBuffer);
        MatrixND<T> DIISExtrapolation(MatrixBuffer& fockMatrices, DIISMatrix& DIISMat);
        void updateWaves(const MatrixND<T>& inv_cholesky);
        void updateSelfConsistentEnergy(DenseVector<T, EDIISBufferSize>& energyBuffer);
    };

    template<class BaseSetType>
    RHFSolver<BaseSetType>::RHFSolver(const Molecular<T>& m, const ElectronConfig& electronConfig_, size_t baseSetSize)
            : molecular(m)
            , electronConfig(electronConfig_)
            , numOccupiedOrbit(electronConfig.getNumOccupiedOrbit())
            , singleHamilton(baseSetSize)
            , overlap(baseSetSize, baseSetSize)
            , baseSet(baseSetSize)
            , selfConsistentEnergy()
            , wave(MatrixND<T>::zeros(baseSetSize, electronConfig.getNumOccupiedOrbit()))
            , eigenSolver(baseSetSize, true) {
        assert(numOccupiedOrbit <= baseSetSize);
    }
    /**
     * Perform self-consistant computation
     * 
     * \return true if converged, false otherwise
     */
    template<class BaseSetType>
    bool RHFSolver<BaseSetType>::compute(const T& criteria, size_t maxIte) {
        assert(criteria.isPositive());

        const size_t baseSetSize = getBaseSetSize();
        formSingleHamilton();
        formOverlapMatrix();

        const MatrixND<T> cholesky = Cholesky(overlap);
        const MatrixND<T> inv_cholesky = cholesky.inv();

        auto densityMatrices = EDIISBuffer(EDIISBufferSize, baseSetSize, baseSetSize, T(0));
        auto sameSpinElectronDensity = MatrixND<T>::zeros(baseSetSize);
        auto fockMatrices = MatrixBuffer(MatrixBufferSize, baseSetSize, baseSetSize, T(0));
        MatrixND<T> fock;
        auto errorMatrices = DIISBuffer(DIISBufferSize - 1, baseSetSize, baseSetSize, T(0));
        DIISMatrix DIISMat = DIISMatrix(DIISBufferSize, DIISBufferSize, -T(1));
        DIISMat[0, 0] = T(0);
        DenseVector<T, EDIISBufferSize> energyBuffer{};

        iteration = 0;
        while (true) {
            auto& abs_error = fock;
            abs_error = abs_elem(*errorMatrices.crbegin());
            const bool nearConverge = abs_error.max() <= T(1E-1); // 1E-1 seleted based on experiments
            const bool doEDIIS = iteration > 0 && (iteration % EDIISBufferSize == 0) && !nearConverge;
            if (doEDIIS)
                EDIISInterpolation(fockMatrices, densityMatrices, energyBuffer);
            else {
                formDensityMatrix(densityMatrices, sameSpinElectronDensity);
                formFockMatrix(fockMatrices, *densityMatrices.crbegin(), sameSpinElectronDensity);
            }

            preDIIS(fockMatrices, errorMatrices, *densityMatrices.crbegin(), inv_cholesky, DIISMat);
            const bool doDIIS = nearConverge && iteration >= DIISBufferSize - 1;
            if (doDIIS)
                fock = DIISExtrapolation(fockMatrices, DIISMat);
            else
                fock = *fockMatrices.crbegin();

            const MatrixND<T> modifiedFock = (inv_cholesky * fock).compute() * inv_cholesky.transpose();
            eigenSolver.compute(modifiedFock);
            eigenSolver.sort();

            updateWaves(inv_cholesky);
            updateSelfConsistentEnergy(energyBuffer);

            if (iteration > 0) {
                const T oldSelfConsistentEnergy = energyBuffer[EDIISBufferSize - 2];
                const T delta = abs((oldSelfConsistentEnergy - selfConsistentEnergy) / oldSelfConsistentEnergy);
                if (delta < criteria)
                    return true;
            }

            if ((++iteration) > maxIte)
                return false;
        }
        return true;
    }

    template<class BaseSetType>
    auto&& RHFSolver<BaseSetType>::getBaseSet(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).baseSet;
    }

    template<class BaseSetType>
    void RHFSolver<BaseSetType>::formSingleHamilton() {
        const size_t baseSetSize = getBaseSetSize();
        for (size_t i = 0; i < baseSetSize; ++i) {
            for (size_t j = i; j < baseSetSize; ++j) {
                T temp = T(0);
                for (size_t k = 0; k < molecular.getAtomCount(); ++k)
                    temp -= BaseSetType::nuclearAttraction(baseSet[i], baseSet[j], molecular.getAtom(k).v())
                            * T(molecular.getAtomicNumber(k));
                singleHamilton[i, j] = BaseSetType::kinetic(baseSet[i], baseSet[j]) + temp;
            }
        }
    }

    template<class BaseSetType>
    void RHFSolver<BaseSetType>::formOverlapMatrix() {
        const size_t baseSetSize = getBaseSetSize();
        for (size_t i = 0; i < baseSetSize; ++i) {
            size_t j = 0;
            for (; j < i; ++j)
                overlap[j, i] = overlap[i, j];

            for (; j < baseSetSize; ++j)
                overlap[j, i] = BaseSetType::overlap(baseSet[i], baseSet[j]);
        }
    }

    template<class BaseSetType>
    void RHFSolver<BaseSetType>::formDensityMatrix(EDIISBuffer& densityMatrices,
                                                   MatrixND<T>& sameSpinElectronDensity) {
        for (size_t i = 0; i < densityMatrices.getLength() - 1; ++i)
            densityMatrices[i].swap(densityMatrices[i + 1]);

        const size_t baseSetSize = getBaseSetSize();
        MatrixND<T>& electronDensity = *densityMatrices.rbegin();
        for (size_t i = 0; i < baseSetSize; ++i) {
            size_t j = 0;
            for (; j < i; ++j) {
                electronDensity[j, i] = electronDensity[i, j];
                sameSpinElectronDensity[j, i] = sameSpinElectronDensity[i, j];
            }

            for (; j < baseSetSize; ++j) {
                T temp1 = T(0);
                T temp2 = T(0);
                for (size_t k = 0; k < numOccupiedOrbit; ++k) {
                    const size_t orbitPos = electronConfig.getOccupiedOrbitPos(k);
                    const auto orbitState = electronConfig.getOrbitState(orbitPos);
                    assert(orbitState != ElectronConfig::NoOccupacy);
                    const bool isSingleOccupacy = orbitState == ElectronConfig::SingleOccupacy;
                    auto orbit = wave.col(k);
                    const T dot = orbit[i] * orbit[j];
                    temp1 += isSingleOccupacy ? dot : T(2) * dot;
                    temp2 += dot;
                }
                electronDensity[j, i] = temp1;
                sameSpinElectronDensity[j, i] = temp2;
            }
        }
    }

    template<class BaseSetType>
    void RHFSolver<BaseSetType>::formFockMatrix(MatrixBuffer& fockMatrices,
                                                const MatrixND<T>& electronDensity,
                                                const MatrixND<T>& sameSpinElectronDensity) {
        const size_t size = getBaseSetSize();
        auto& fock = *fockMatrices.begin();
        for (size_t p = 0; p < size; ++p) {
            for (size_t q = 0; q < size; ++q) {
                T temp = T(0);
                for (size_t r = 0; r < size; ++r) {
                    for (size_t s = 0; s < size; ++s) {
                        const T coulomb = BaseSetType::electronRepulsion(baseSet[p], baseSet[r], baseSet[q], baseSet[s]);
                        const T exchange = BaseSetType::electronRepulsion(baseSet[p], baseSet[r], baseSet[s], baseSet[q]);
                        temp += electronDensity[s, r] * coulomb - sameSpinElectronDensity[s, r] * exchange;
                    }
                }
                fock[q, p] = temp;
            }
        }
        fock += singleHamilton;
        for (size_t i = 0; i < fockMatrices.getLength() - 1; ++i)
            fockMatrices[i].swap(fockMatrices[i + 1]);
    }

    template<class BaseSetType>
    void RHFSolver<BaseSetType>::preDIIS(const MatrixBuffer& fockMatrices,
                                         DIISBuffer& errorMatrices,
                                         const MatrixND<T>& electronDensity,
                                         const MatrixND<T>& inv_cholesky,
                                         DIISMatrix& DIISMat) {
        /* Insert next error matrix */ {
            const MatrixND<T> term1 = (*fockMatrices.crbegin() * electronDensity).compute() * overlap;
            const MatrixND<T> term2 = (overlap * electronDensity).compute() * (*fockMatrices.crbegin());
            const MatrixND<T> temp = term1 - term2;
            errorMatrices[0] = (inv_cholesky * temp).compute() * inv_cholesky.transpose();
            for (size_t i = 0; i < errorMatrices.getLength() - 1; ++i)
                errorMatrices[i].swap(errorMatrices[i + 1]);
        }
        /* Construct equation */ {
            for (size_t i = 1; i < DIISMat.getRow(); ++i) {
                for (size_t j = i; j < DIISMat.getRow(); ++j) {
                    T temp = (errorMatrices[i - 1] * errorMatrices[j - 1]).trace();
                    DIISMat[i, j] = temp;
                    DIISMat[j, i] = temp;
                }
            }
        }
    }

    template<class BaseSetType>
    void RHFSolver<BaseSetType>::EDIISInterpolation(MatrixBuffer& fockMatrices,
                                                    EDIISBuffer& densityMatrices,
                                                    const DenseVector<T, EDIISBufferSize>& energyBuffer) {
        constexpr size_t problemDim = EDIISBuffer::getLength();
        auto G = DenseSymmMatrix<T>(problemDim, T(0));
        MatrixND<T> deltaFock;
        MatrixND<T> deltaDensity;
        constexpr size_t offset = MatrixBuffer::getLength() - problemDim;
        for (size_t r = 0; r < problemDim; ++r) {
            for (size_t c = r + 1; c < problemDim; ++c) {
                deltaFock = fockMatrices[offset + r] - fockMatrices[offset + c];
                deltaDensity = densityMatrices[c] - densityMatrices[r];
                G[r, c] = (deltaFock * deltaDensity).trace();
            }
        }

        auto equalityConstraint = DenseMatrix<T, MatrixOption::Row, 1, Dynamic>(1, problemDim + 1, T(1));
        auto inequalityConstraint = DenseMatrix<T, MatrixOption::Row>(problemDim, problemDim + 1, T(0));
        auto block = inequalityConstraint.leftCols(problemDim);
        block.toIdentity();
        auto initial = VectorND<T>(problemDim, T(0));
        QuadraticProgramming<T> QP(G, energyBuffer, equalityConstraint, inequalityConstraint, initial);
        QP.compute_nonconvex(T(1E-10));

        auto& newFock = deltaFock;
        auto& newDensity = deltaDensity;
        newFock = T(0);
        newDensity = T(0);
        for (size_t i = 0; i < problemDim; ++i) {
            newFock += fockMatrices[offset + i] * QP.getSolution()[i];
            newDensity += densityMatrices[i] * QP.getSolution()[i];
        }

        for (size_t i = 0; i < fockMatrices.getLength() - 1; ++i)
            fockMatrices[i].swap(fockMatrices[i + 1]);
        for (size_t i = 0; i < densityMatrices.getLength() - 1; ++i)
            densityMatrices[i].swap(densityMatrices[i + 1]);
        *fockMatrices.rbegin() = newFock;
        *densityMatrices.rbegin() = newDensity;
    }

    template<class BaseSetType>
    auto RHFSolver<BaseSetType>::DIISExtrapolation(MatrixBuffer& fockMatrices, DIISMatrix& DIISMat) -> MatrixND<T> {
        DenseVector<T, DIISBufferSize> x{};
        /* Solve linear equation */ {
            DenseVector<T, DIISBufferSize> b = DenseVector<T, DIISBufferSize>(DIISBufferSize, T(0));
            b[0] = -T(1);
            const DIISMatrix inv_A = DIISMat.inv();
            x = inv_A * b;
        }

        auto extrapolate_fock = MatrixND<T>::zeros(getBaseSetSize());
        constexpr size_t offset = MatrixBuffer::getLength() - DIISBuffer::getLength();
        for (size_t i = 1; i < x.getLength(); ++i)
            extrapolate_fock += fockMatrices[offset + i - 1] * x[i];
        return extrapolate_fock;
    }

    template<class BaseSetType>
    void RHFSolver<BaseSetType>::updateWaves(const MatrixND<T>& inv_cholesky) {
        const auto& eigenvectors = eigenSolver.getRawEigenvectors();
        for (size_t i = 0; i < numOccupiedOrbit; ++i) {
            auto eigenState = wave.col(i);
            const size_t orbitPos = electronConfig.getOccupiedOrbitPos(i);
            eigenState = inv_cholesky.transpose() * eigenvectors.col(orbitPos);
        }
    }

    template<class BaseSetType>
    void RHFSolver<BaseSetType>::updateSelfConsistentEnergy(DenseVector<T, EDIISBufferSize>& energyBuffer) {
        for (size_t i = 0; i < energyBuffer.getLength() - 1; ++i)
            energyBuffer[i].swap(energyBuffer[i + 1]);

        const auto& eigenvalues = eigenSolver.getEigenvalues();
        selfConsistentEnergy = T(0);
        for (size_t i = 0; i < wave.getCol(); ++i) {
            const size_t orbitPos = electronConfig.getOccupiedOrbitPos(i);
            T temp = eigenvalues[orbitPos].real();
            auto orbit = wave.col(i);
            temp += (orbit.transpose() * singleHamilton).compute().row(0) * orbit;
            const auto orbitState = electronConfig.getOrbitState(orbitPos);
            assert(orbitState != ElectronConfig::NoOccupacy);
            const bool isSingleOccupacy = orbitState == ElectronConfig::SingleOccupacy;
            selfConsistentEnergy += isSingleOccupacy ? temp : (T(2) * temp);
        }
        selfConsistentEnergy *= T(0.5);
        energyBuffer[EDIISBufferSize - 1] = selfConsistentEnergy;
    }
}
