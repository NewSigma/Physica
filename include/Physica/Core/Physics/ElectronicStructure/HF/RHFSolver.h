/*
 * Copyright 2021-2022 WeiBo He.
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
#include "Physica/Core/Physics/ElectronicStructure/ElectronConfig.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/Cholesky.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/InverseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Transpose.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/EigenSolver.h"
#include "Physica/Core/Math/Optimization/QuadraticProgramming/QuadraticProgramming.h"

namespace Physica::Core::Physics {
    namespace Internal {
        template<class T> class Traits;
    }
    /**
     * Reference:
     * [1] Jos Thijssen. Computational Physics[M].London: Cambridge university press, 2013:43-88
     * [2] Larsen A, Poulsen R S. Applied Hartree-Fock methods.
     * [3] Kudin K N, Scuseria G E, Cances E. A black-box self-consistent field convergence algorithm: One step closer[J]. Journal of Chemical Physics, 2002, 116(19):8255-8261.
     */
    template<class BaseSetType>
    class RHFSolver {
        using ScalarType = typename Internal::Traits<BaseSetType>::ScalarType;
        using MatrixType = DenseMatrix<ScalarType, DenseMatrixOption::Column | DenseMatrixOption::Vector>;

        constexpr static size_t EDIISBufferSize = 3; //Refer EDIIS from [3]
        constexpr static size_t DIISBufferSize = 3; //Refer DIIS from [2]
        constexpr static size_t MatrixBufferSize = std::max(EDIISBufferSize, DIISBufferSize);
        static_assert(DIISBufferSize >= 3, "DIISBufferSize less than three makes no sence");
        using EDIISBuffer = Utils::Array<MatrixType, EDIISBufferSize>;
        using DIISBuffer = Utils::Array<MatrixType, DIISBufferSize - 1>;
        using MatrixBuffer = Utils::Array<MatrixType, MatrixBufferSize>;
        using DIISMatrix = DenseMatrix<ScalarType, DenseMatrixOption::Column | DenseMatrixOption::Element, DIISBufferSize, DIISBufferSize>;
    public:
        using WaveType = MatrixType;
    private:
        const Molecular<ScalarType>& molecular;
        ElectronConfig electronConfig;
        size_t numOccupiedOrbit;
        DenseSymmMatrix<ScalarType, Dynamic> singleHamilton;
        MatrixType overlap;
        Utils::Array<BaseSetType> baseSet;
        ScalarType selfConsistentEnergy;
        MatrixType wave;
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
        [[nodiscard]] const MatrixType& getWave() const noexcept { return wave; }
        [[nodiscard]] size_t getIteration() const noexcept { return iteration; }
        /* Setters */
        void setInitialWave(const MatrixType& initialWave) { wave = initialWave; }
    private:
        void formSingleHamilton();
        void formOverlapMatrix();
        void formDensityMatrix(EDIISBuffer& densityMatrices,
                               MatrixType& sameSpinElectronDensity);
        void formFockMatrix(MatrixBuffer& fockMatrices,
                            const MatrixType& electronDensity,
                            const MatrixType& sameSpinElectronDensity);
        void preDIIS(const MatrixBuffer& fockMatrices,
                     DIISBuffer& errorMatrices,
                     const MatrixType& electronDensity,
                     const MatrixType& inv_cholesky,
                     DIISMatrix& DIISMat);
        void EDIISInterpolation(MatrixBuffer& fockMatrices,
                                EDIISBuffer& densityMatrices,
                                const Vector<ScalarType, EDIISBufferSize>& energyBuffer);
        MatrixType DIISExtrapolation(MatrixBuffer& fockMatrices, DIISMatrix& DIISMat);
        void updateWaves(const MatrixType& inv_cholesky);
        [[nodiscard]] ScalarType updateSelfConsistentEnergy(Vector<ScalarType, EDIISBufferSize>& energyBuffer);
    };

    template<class BaseSetType>
    RHFSolver<BaseSetType>::RHFSolver(const Molecular<ScalarType>& m, const ElectronConfig& electronConfig_, size_t baseSetSize)
            : molecular(m)
            , electronConfig(electronConfig_)
            , numOccupiedOrbit(electronConfig.getNumOccupiedOrbit())
            , singleHamilton(baseSetSize)
            , overlap(baseSetSize, baseSetSize)
            , baseSet(baseSetSize)
            , selfConsistentEnergy()
            , wave(MatrixType::Zeros(baseSetSize, electronConfig.getNumOccupiedOrbit()))
            , eigenSolver(baseSetSize)
            , iteration(0) {
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

        auto densityMatrices = EDIISBuffer(EDIISBufferSize, MatrixType::Zeros(baseSetSize, baseSetSize));
        MatrixType sameSpinElectronDensity = MatrixType::Zeros(baseSetSize);
        auto fockMatrices = MatrixBuffer(MatrixBufferSize, MatrixType::Zeros(baseSetSize, baseSetSize));
        MatrixType fock;
        auto errorMatrices = DIISBuffer(DIISBufferSize - 1, MatrixType::Zeros(baseSetSize, baseSetSize));
        DIISMatrix DIISMat = DIISMatrix(DIISBufferSize, DIISBufferSize, -ScalarType::One());
        DIISMat(0, 0) = ScalarType::Zero();
        Vector<ScalarType, EDIISBufferSize> energyBuffer{};

        iteration = 0;
        do {
            MatrixType& abs_error = fock;
            abs_error = abs(*errorMatrices.crbegin());
            const bool nearConverge = abs_error.max() <= ScalarType(1E-1); //1E-1 seleted by experiment
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

            const MatrixType modifiedFock = (inv_cholesky * fock).compute() * inv_cholesky.transpose();
            eigenSolver.compute(modifiedFock, true);
            eigenSolver.sort();

            updateWaves(inv_cholesky);

            const ScalarType delta = updateSelfConsistentEnergy(energyBuffer);
            if (delta < criteria)
                return true;

            if ((++iteration) > maxIte)
                return false;
        } while(true);
        return true;
    }

    template<class BaseSetType>
    void RHFSolver<BaseSetType>::formSingleHamilton() {
        const size_t baseSetSize = getBaseSetSize();
        for (size_t i = 0; i < baseSetSize; ++i) {
            for (size_t j = i; j < baseSetSize; ++j) {
                ScalarType temp = ScalarType::Zero();
                for (size_t k = 0; k < molecular.getAtomCount(); ++k)
                    temp -= BaseSetType::nuclearAttraction(baseSet[i], baseSet[j], molecular.getAtom(k).v())
                            * ScalarType(molecular.getAtomicNumber(k));
                singleHamilton(i, j) = BaseSetType::kinetic(baseSet[i], baseSet[j]) + temp;
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
    void RHFSolver<BaseSetType>::formDensityMatrix(EDIISBuffer& densityMatrices,
                                                   MatrixType& sameSpinElectronDensity) {
        for (size_t i = 0; i < densityMatrices.getLength() - 1; ++i)
            swap(densityMatrices[i], densityMatrices[i + 1]);

        const size_t baseSetSize = getBaseSetSize();
        MatrixType& electronDensity = *densityMatrices.rbegin();
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
                    auto orbit = wave.col(k);
                    const ScalarType dot = orbit[i] * orbit[j];
                    temp1 += isSingleOccupacy ? dot : ScalarType::Two() * dot;
                    temp2 += dot;
                }
                electronDensity(j, i) = temp1;
                sameSpinElectronDensity(j, i) = temp2;
            }
        }
    }

    template<class BaseSetType>
    void RHFSolver<BaseSetType>::formFockMatrix(MatrixBuffer& fockMatrices,
                                                const MatrixType& electronDensity,
                                                const MatrixType& sameSpinElectronDensity) {
        const size_t size = getBaseSetSize();
        auto& fock = *fockMatrices.begin();
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
        fock += singleHamilton;
        for (size_t i = 0; i < fockMatrices.getLength() - 1; ++i)
            swap(fockMatrices[i], fockMatrices[i + 1]);
    }

    template<class BaseSetType>
    void RHFSolver<BaseSetType>::preDIIS(const MatrixBuffer& fockMatrices,
                                         DIISBuffer& errorMatrices,
                                         const MatrixType& electronDensity,
                                         const MatrixType& inv_cholesky,
                                         DIISMatrix& DIISMat) {
        /* Insert next error matrix */ {
            const MatrixType term1 = (*fockMatrices.crbegin() * electronDensity).compute() * overlap;
            const MatrixType term2 = (overlap * electronDensity).compute() * (*fockMatrices.crbegin());
            const MatrixType temp = term1 - term2;
            errorMatrices[0] = (inv_cholesky * temp).compute() * inv_cholesky.transpose();
            for (size_t i = 0; i < errorMatrices.getLength() - 1; ++i)
                swap(errorMatrices[i], errorMatrices[i + 1]);
        }
        /* Construct equation */ {
            for (size_t i = 1; i < DIISMat.getRow(); ++i) {
                for (size_t j = i; j < DIISMat.getRow(); ++j) {
                    ScalarType temp = (errorMatrices[i - 1] * errorMatrices[j - 1]).trace();
                    DIISMat(i, j) = temp;
                    DIISMat(j, i) = temp;
                }
            }
        }
    }

    template<class BaseSetType>
    void RHFSolver<BaseSetType>::EDIISInterpolation(MatrixBuffer& fockMatrices,
                                                    EDIISBuffer& densityMatrices,
                                                    const Vector<ScalarType, EDIISBufferSize>& energyBuffer) {
        constexpr size_t problemDim = EDIISBuffer::getLength();
        auto G = DenseSymmMatrix<ScalarType>(problemDim, ScalarType::Zero());
        MatrixType deltaFock;
        MatrixType deltaDensity;
        constexpr size_t offset = MatrixBuffer::getLength() - problemDim;
        for (size_t r = 0; r < problemDim; ++r) {
            for (size_t c = r + 1; c < problemDim; ++c) {
                deltaFock = fockMatrices[offset + r] - fockMatrices[offset + c];
                deltaDensity = densityMatrices[c] - densityMatrices[r];
                G(r, c) = (deltaFock * deltaDensity).trace();
            }
        }

        auto equalityConstraint = DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Element, 1, Dynamic>(1, problemDim + 1, ScalarType::One());
        auto inequalityConstraint = DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Vector>(problemDim, problemDim + 1, ScalarType::Zero());
        auto block = inequalityConstraint.leftCols(problemDim);
        block.toUnitMatrix();
        auto initial = Vector<ScalarType>(problemDim, ScalarType::Zero());
        QuadraticProgramming<ScalarType> QP(G, energyBuffer, equalityConstraint, inequalityConstraint, initial);
        QP.compute_nonconvex(ScalarType(1E-10));

        auto& newFock = deltaFock;
        auto& newDensity = deltaDensity;
        newFock = ScalarType::Zero();
        newDensity = ScalarType::Zero();
        for (size_t i = 0; i < problemDim; ++i) {
            newFock += fockMatrices[offset + i] * QP.getSolution()[i];
            newDensity += densityMatrices[i] * QP.getSolution()[i];
        }

        for (size_t i = 0; i < fockMatrices.getLength() - 1; ++i)
            swap(fockMatrices[i], fockMatrices[i + 1]);
        for (size_t i = 0; i < densityMatrices.getLength() - 1; ++i)
            swap(densityMatrices[i], densityMatrices[i + 1]);
        *fockMatrices.rbegin() = newFock;
        *densityMatrices.rbegin() = newDensity;
    }

    template<class BaseSetType>
    typename RHFSolver<BaseSetType>::MatrixType RHFSolver<BaseSetType>::DIISExtrapolation(MatrixBuffer& fockMatrices,
                                                                                          DIISMatrix& DIISMat) {
        Vector<ScalarType, DIISBufferSize> x{};
        /* Solve linear equation */ {
            Vector<ScalarType, DIISBufferSize> b = Vector<ScalarType, DIISBufferSize>(DIISBufferSize, ScalarType::Zero());
            b[0] = -ScalarType::One();
            const DIISMatrix inv_A = DIISMat.inverse();
            x = inv_A * b;
        }

        MatrixType extrapolate_fock = MatrixType::Zeros(getBaseSetSize());
        constexpr size_t offset = MatrixBuffer::getLength() - DIISBuffer::getLength();
        for (size_t i = 1; i < x.getLength(); ++i)
            extrapolate_fock += fockMatrices[offset + i - 1] * x[i];
        return extrapolate_fock;
    }

    template<class BaseSetType>
    void RHFSolver<BaseSetType>::updateWaves(const MatrixType& inv_cholesky) {
        const auto& eigenvectors = eigenSolver.getRawEigenvectors();
        for (size_t i = 0; i < numOccupiedOrbit; ++i) {
            auto eigenState = wave.col(i);
            const size_t orbitPos = electronConfig.getOccupiedOrbitPos(i);
            eigenState = inv_cholesky.transpose() * eigenvectors.col(orbitPos);
        }
    }

    template<class BaseSetType>
    typename RHFSolver<BaseSetType>::ScalarType RHFSolver<BaseSetType>::updateSelfConsistentEnergy(
            Vector<ScalarType, EDIISBufferSize>& energyBuffer) {
        for (size_t i = 0; i < energyBuffer.getLength() - 1; ++i)
            swap(energyBuffer[i], energyBuffer[i + 1]);

        const auto& eigenvalues = eigenSolver.getEigenvalues();
        selfConsistentEnergy = ScalarType::Zero();
        for (size_t i = 0; i < wave.getColumn(); ++i) {
            const size_t orbitPos = electronConfig.getOccupiedOrbitPos(i);
            ScalarType temp = eigenvalues[orbitPos].getReal();
            auto orbit = wave.col(i);
            temp += (orbit.asMatrix().transpose() * singleHamilton).compute().row(0) * orbit;
            const auto orbitState = electronConfig.getOrbitState(orbitPos);
            assert(orbitState != ElectronConfig::NoOccupacy);
            const bool isSingleOccupacy = orbitState == ElectronConfig::SingleOccupacy;
            selfConsistentEnergy += isSingleOccupacy ? temp : (ScalarType::Two() * temp);
        }
        selfConsistentEnergy *= ScalarType(0.5);
        auto ite = energyBuffer.rbegin();
        *ite = selfConsistentEnergy;
        const ScalarType oldSelfConsistentEnergy = *(++ite);
        const ScalarType delta = abs((oldSelfConsistentEnergy - selfConsistentEnergy) / oldSelfConsistentEnergy);
        return delta;
    }
}
