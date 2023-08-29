/*
 * Copyright 2023 WeiBo He.
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

#include <iostream>
#include "Physica/Core/Physics/MD/MDCell.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Grid/RSpaceGrid.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/EigenSolver.h"
#include "Physica/Core/Math/Calculus/Interpolation.h"
#include "Physica/Core/Math/Calculus/PDE/FEM/Element/CuboidLinear.h"
#include "Physica/Core/Math/Statistics/NumCharacter.h"
#include "Physica/Core/Math/Transform/FFT.h"
#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"

namespace Physica::Core {
    /**
     * References:
     * [1] Dario Alfè PHON: A program to calculate phonons using the small displacement method [J]. Computer Physics Communications, 2009, 180(12), 2622-2633
     */
    template<class ScalarType, class PosScalarType>
    class FrozenPhonon {
        using This = FrozenPhonon<ScalarType, PosScalarType>;
        using ComplexType = ComplexScalar<ScalarType>;
        using Index3D = Utils::Array<size_t, 3>;
        using Vector3D = Vector<ScalarType, 3>;
        using FFT3D = FFT<ScalarType, 3>;
    public:
        using MDCellType = MDCell<ScalarType, PosScalarType>;
        using PositionMatrix = typename MDCellType::PositionMatrix;
        using MatrixType = DenseMatrix<ComplexType>;
        using MatrixGrid = GridStorage<MatrixType>;
        using EigenSolverType = EigenSolver<MatrixType>;
        using QPointGrid = GridStorage<EigenSolverType>;
        constexpr static size_t Dim = Internal::Traits<MDCellType>::Dim;
    private:
        MDCellType unitCell;
        Index3D superSize;
    public:
        FrozenPhonon(MDCellType unitCell_, Index3D superSize_);
        FrozenPhonon(const FrozenPhonon&) = default;
        FrozenPhonon(FrozenPhonon&&) noexcept = default;
        ~FrozenPhonon() = default;
        /* Operators */
        FrozenPhonon& operator=(FrozenPhonon obj) noexcept;
        /* Operations */
        template<class ForceModel>
        [[nodiscard]] MatrixGrid makeForceConstants(const ForceModel& model, ScalarType displace, ScalarType translationPrec, size_t maxIteration) const;
        void toDynamicMatrix(MatrixType& forceConstant) const;
        void toDynamicMatrix(MatrixGrid& forceConstants) const;
        void swap(FrozenPhonon& obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumUnitCellAtom() const noexcept { return unitCell.getNumParticle(); }
        [[nodiscard]] size_t getUnitCellDOF() const noexcept { return Dim * getNumUnitCellAtom(); }
        [[nodiscard]] size_t getNumSuperCellAtom() const noexcept { return getNumUnitCellAtom() * getNumCell(); }
        [[nodiscard]] size_t getSuperCellDOF() const noexcept { return getUnitCellDOF() * getNumCell(); }        [[nodiscard]] Index3D getSuperSize() const noexcept { return superSize; }
        [[nodiscard]] Index3D getForceConstantsGridSize() const noexcept { return FFT3D::rSizeToKSize(superSize); }
        [[nodiscard]] size_t getNumCell() const noexcept { return superSize[0] * superSize[1] * superSize[2]; }
        [[nodiscard]] Vector<ScalarType> makeFreq(const EigenSolverType& eigen) const;
        [[nodiscard]] inline Vector<ScalarType> makeFreq(const QPointGrid& qPoints, Index3D qIndex) const;
        [[nodiscard]] DenseMatrix<ScalarType> makeEigenVectors(const EigenSolverType& eigen) const;
        [[nodiscard]] inline DenseMatrix<ScalarType> makeEigenVectors(const QPointGrid& qPoints, Index3D qIndex) const;
        /* Static members */
        [[nodiscard]] static EigenSolverType diagonalize(const MatrixType& dynamicMatrix);
        [[nodiscard]] static QPointGrid diagonalize(const MatrixGrid& dynamicMatrixes);
    private:
        ScalarType removeDriftForce(Vector<ScalarType>& force) const;
    };

    template<class ScalarType, class PosScalarType>
    FrozenPhonon<ScalarType, PosScalarType>::FrozenPhonon(MDCellType unitCell_, Index3D superSize_)
            : unitCell(std::move(unitCell_))
            , superSize(superSize_) {}

    template<class ScalarType, class PosScalarType>
    FrozenPhonon<ScalarType, PosScalarType>& FrozenPhonon<ScalarType, PosScalarType>::operator=(FrozenPhonon obj) noexcept {
        swap(obj);
        return *this;
    }
    /**
     * Use iteration method to apply translational invariance and while keep force constant matrix symmetric as introduced in [1].
     */
    template<class ScalarType, class PosScalarType>
    template<class ForceModel>
    typename FrozenPhonon<ScalarType, PosScalarType>::MatrixGrid
    FrozenPhonon<ScalarType, PosScalarType>::makeForceConstants(
            const ForceModel& model,
            ScalarType displace,
            ScalarType translationPrec,
            size_t maxIteration) const {
        const size_t unitCellDOF = getUnitCellDOF();
        const ScalarType factor = -reciprocal(displace);
        const MDCellType superCell = unitCell.template makeSuperCell<ExtendCellOption::CellMajor>(superSize);
        
        MatrixGrid result(getForceConstantsGridSize(), unitCellDOF, unitCellDOF, ScalarType(0));
        auto& fcMatrixes = result.asArray();
        FFT3D fft(superSize, {1, 1, 1}, PlanFlag::Estimate);
        auto rSpace = fft.getRSpace().flatten();
        auto kSpace = fft.getKSpace().flatten();
        /* Make symmetrize force constants matrix */ {
            PositionMatrix pos = superCell.getPos();
            for (size_t row = 0; row < unitCellDOF; ++row) {
                ScalarType& toDisplace = pos(row / Dim, row % Dim);
                const ScalarType copy = toDisplace;
                toDisplace += displace;
                Vector<ScalarType> forceConst =
                        model.template force<SequentialExecutor>(MDCellType(superCell.getLattice(), pos, superCell.getMassVec())) * factor;
                toDisplace = copy;

                for (size_t col = 0; col < unitCellDOF; ++col) {
                    const size_t shift = unitCellDOF;
                    for (size_t cell = 0; cell < getNumCell(); ++cell)
                        rSpace[cell] = forceConst[col + cell * shift];
                    fft.transform();

                    for (size_t i = 0; i < fcMatrixes.getLength(); ++i) {
                        auto& fcMatrix = fcMatrixes[i];
                        const ComplexType elem = kSpace[i] * ScalarType(0.5);
                        fcMatrix(row, col) += elem;
                        fcMatrix(col, row) += elem.conjugate();
                    }
                }
            }
        }
        /* Apply translational invariance */ {
            Vector<ScalarType> forceConst(getSuperCellDOF());
            MatrixType temp;
            ScalarType averageDrift = std::numeric_limits<ScalarType>::max();
            size_t iteration = 0;
            while (averageDrift > translationPrec) {
                if (iteration == maxIteration) [[unlikely]]
                    throw BadConvergenceException("[Error]: Failed to apply translational invariance in the given steps");
                iteration += 1;

                for (size_t row = 0; row < unitCellDOF; ++row) {
                    for (size_t col = 0; col < unitCellDOF; ++col) {
                        for (size_t i = 0; i < fcMatrixes.getLength(); ++i) {
                            auto& fcMatrix = fcMatrixes[i];
                            kSpace[i] = fcMatrix(row, col);
                        }
                        fft.invTransform();

                        const size_t shift = unitCellDOF;
                        for (size_t cell = 0; cell < getNumCell(); ++cell)
                            forceConst[col + cell * shift] = rSpace[cell];
                    }
                    toNextMean(averageDrift, row, removeDriftForce(forceConst));

                    for (size_t col = 0; col < unitCellDOF; ++col) {
                        const size_t shift = unitCellDOF;
                        for (size_t cell = 0; cell < getNumCell(); ++cell)
                            rSpace[cell] = forceConst[col + cell * shift];
                        fft.transform();

                        for (size_t i = 0; i < fcMatrixes.getLength(); ++i) {
                            auto& fcMatrix = fcMatrixes[i];
                            fcMatrix(row, col) = kSpace[i];
                        }
                    }
                }

                for (size_t i = 0; i < fcMatrixes.getLength(); ++i) {
                    auto& fcMatrix = fcMatrixes[i];
                    temp = (fcMatrix + fcMatrix.transpose().conjugate()) * ScalarType(0.5);
                    fcMatrix = temp;
                }
            }
        }
        return result;
    }

    template<class ScalarType, class PosScalarType>
    void FrozenPhonon<ScalarType, PosScalarType>::toDynamicMatrix(MatrixType& forceConstant) const {
        const size_t unitCellDOF = getUnitCellDOF();
        for (size_t row = 0; row < unitCellDOF; ++row) {
            const size_t atom1 = row / Dim;
            const ScalarType mass1 = unitCell.getMass(atom1);
            for (size_t col = 0; col < unitCellDOF; ++col) {
                const size_t atom2 = col / Dim;
                const ScalarType mass2 = unitCell.getMass(atom2);
                const ScalarType repMass = reciprocal(sqrt(mass1 * mass2));
                forceConstant(row, col) *= repMass;
            }
        }
    }

    template<class ScalarType, class PosScalarType>
    void FrozenPhonon<ScalarType, PosScalarType>::toDynamicMatrix(MatrixGrid& forceConstants) const {
        const size_t unitCellDOF = getUnitCellDOF();
        auto& fcMatrixes = forceConstants.flatten();
        for (size_t row = 0; row < unitCellDOF; ++row) {
            const size_t atom1 = row / Dim;
            const ScalarType mass1 = unitCell.getMass(atom1);
            for (size_t col = 0; col < unitCellDOF; ++col) {
                const size_t atom2 = col / Dim;
                const ScalarType mass2 = unitCell.getMass(atom2);
                const ScalarType repMass = reciprocal(sqrt(mass1 * mass2));

                for (size_t i = 0; i < fcMatrixes.getLength(); ++i) {
                    auto& fcMatrix = fcMatrixes[i];
                    fcMatrix(row, col) *= repMass;
                }
            }
        }
    }

    template<class ScalarType, class PosScalarType>
    void FrozenPhonon<ScalarType, PosScalarType>::swap(This& obj) noexcept {
        unitCell.swap(obj.unitCell);
        superSize.swap(obj.superSize);
    }

    template<class ScalarType, class PosScalarType>
    Vector<ScalarType> FrozenPhonon<ScalarType, PosScalarType>::makeFreq(const EigenSolverType& eigen) const {
        Vector<ScalarType> result(getUnitCellDOF());
        const Vector<ScalarType> eigenvalues = toRealVector(eigen.getEigenvalues());
        for (size_t i = 0; i < getUnitCellDOF(); ++i)
            result[i] = eigenvalues[i].isNegative() ? -1 : 1;
        result = hadamard(result, sqrt(abs(eigenvalues)));
        result *= ScalarType(1 / (2 * M_PI));
        return result;
    }

    template<class ScalarType, class PosScalarType>
    inline Vector<ScalarType> FrozenPhonon<ScalarType, PosScalarType>::makeFreq(const QPointGrid& qPoints, Index3D qIndex) const {
        return makeFreq(qPoints(qIndex));
    }

    template<class ScalarType, class PosScalarType>
    DenseMatrix<ScalarType> FrozenPhonon<ScalarType, PosScalarType>::makeEigenVectors(const EigenSolverType& eigen) const {
        DenseMatrix<ScalarType> result(getUnitCellDOF(), getUnitCellDOF());
        for (size_t i = 0; i < result.getColumn(); ++i) {
            const auto fromCol = eigen.getRawEigenvectors().col(i);
            auto toCol = result.col(i);
            for (size_t j = 0; j < result.getRow(); ++j) {
                const ScalarType repSqrtMass = reciprocal(sqrt(unitCell.getMass(j / Dim)));
                toCol[j] = fromCol[j].getReal() * repSqrtMass;
            }
            toCol.toUnit(); // Optimize: repSqrtMass matrix is diagonal
        }
        return result;
    }

    template<class ScalarType, class PosScalarType>
    inline DenseMatrix<ScalarType> FrozenPhonon<ScalarType, PosScalarType>::makeEigenVectors(const QPointGrid& qPoints, Index3D qIndex) const {
        return makeEigenVectors(qPoints(qIndex));
    }

    template<class ScalarType, class PosScalarType>
    typename FrozenPhonon<ScalarType, PosScalarType>::EigenSolverType
    FrozenPhonon<ScalarType, PosScalarType>::diagonalize(const MatrixType& dynamicMatrix) {
        const size_t unitCellDOF = dynamicMatrix.getRow();
        auto eigen = EigenSolverType(unitCellDOF);
        eigen.compute(dynamicMatrix, true);
        eigen.sort();
        return eigen;
    }

    template<class ScalarType, class PosScalarType>
    typename FrozenPhonon<ScalarType, PosScalarType>::QPointGrid FrozenPhonon<ScalarType, PosScalarType>::diagonalize(
            const MatrixGrid& dynamicMatrixes) {
        const auto& matrixes = dynamicMatrixes.flatten();
        const size_t unitCellDOF = matrixes[0].getRow();
        QPointGrid qPoints(dynamicMatrixes.getDim(), unitCellDOF);
        for (size_t i = 0; i < matrixes.getLength(); ++i) {
            auto& eigen = qPoints.flatten()[i];
            eigen.compute(matrixes[i], true);
            eigen.sort();
        }
        return qPoints;
    }

    template<class ScalarType, class PosScalarType>
    ScalarType FrozenPhonon<ScalarType, PosScalarType>::removeDriftForce(Vector<ScalarType>& force) const {
        const size_t superCellDOF = getSuperCellDOF();
        assert(force.getLength() == superCellDOF);
        const size_t numSuperAtom = superCellDOF / Dim;
        const ScalarType factor = reciprocal(ScalarType(numSuperAtom));
        ScalarType averageDrift = 0;
        for (size_t dim = 0; dim < Dim; ++dim) {
            ScalarType drift = 0;
            for (size_t i = dim; i < superCellDOF; i += Dim)
                drift += force[i];
            drift *= factor;
            for (size_t i = dim; i < superCellDOF; i += Dim)
                force[i] -= drift;
            averageDrift += abs(drift);
        }
        averageDrift /= ScalarType(Dim);
        return averageDrift;
    }
}
