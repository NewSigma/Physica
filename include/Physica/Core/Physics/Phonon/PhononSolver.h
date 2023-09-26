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

#include "Physica/Core/Physics/MD/MDCell.h"

namespace Physica::Core {
    template<class ScalarType, class PosScalarType>
    class PhononSolver {
        using This = PhononSolver<ScalarType, PosScalarType>;
    public:
        using ComplexType = ComplexScalar<ScalarType>;
        using Index3D = Utils::Array<size_t, 3>;
        using Vector3D = Vector<ScalarType, 3>;
        using FFT3D = FFT<ScalarType, 3>;
        using MDCellType = MDCell<ScalarType, PosScalarType>;
        using MatrixType = DenseMatrix<ComplexType>;
        using PositionMatrix = typename MDCellType::PositionMatrix;
        using MatrixGrid = GridStorage<MatrixType>;
        using EigenSolverType = EigenSolver<MatrixType>;
        using QPointGrid = GridStorage<EigenSolverType>;
        constexpr static unsigned int Dim = Internal::Traits<MDCellType>::Dim;
    private:
        MDCellType unitCell;
        Index3D superSize;
    public:
        PhononSolver(MDCellType unitCell_, Index3D superSize_);
        PhononSolver(const PhononSolver&) = default;
        PhononSolver(PhononSolver&&) noexcept = default;
        ~PhononSolver() = default;
        /* Operators */
        PhononSolver& operator=(PhononSolver obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] MatrixType interpolatePoint(Vector3D qPoint, const MatrixGrid& forceConstants) const;
        void toDynamicMatrix(MatrixType& forceConstant) const;
        void toDynamicMatrix(MatrixGrid& forceConstants) const;
        [[nodiscard]] Vector<ScalarType> makeFreq(const EigenSolverType& eigen) const;
        [[nodiscard]] inline Vector<ScalarType> makeFreq(const QPointGrid& qPoints, Index3D qIndex) const;
        [[nodiscard]] DenseMatrix<ScalarType> makeEigenVectors(const EigenSolverType& eigen) const;
        [[nodiscard]] inline DenseMatrix<ScalarType> makeEigenVectors(const QPointGrid& qPoints, Index3D qIndex) const;
        void swap(PhononSolver& obj) noexcept;
        /* Getters */
        [[nodiscard]] const MDCellType& getUnitCell() const noexcept { return unitCell; }
        [[nodiscard]] size_t getNumUnitCellAtom() const noexcept { return unitCell.getNumParticle(); }
        [[nodiscard]] size_t getUnitCellDOF() const noexcept { return Dim * getNumUnitCellAtom(); }
        [[nodiscard]] size_t getNumBand() const noexcept { return getUnitCellDOF(); }
        [[nodiscard]] size_t getNumSuperCellAtom() const noexcept { return getNumUnitCellAtom() * getNumCell(); }
        [[nodiscard]] size_t getSuperCellDOF() const noexcept { return getUnitCellDOF() * getNumCell(); }
        [[nodiscard]] Index3D getSuperSize() const noexcept { return superSize; }
        [[nodiscard]] Index3D getForceConstantsGridSize() const noexcept { return FFT3D::rSizeToKSize(superSize); }
        [[nodiscard]] size_t getNumCell() const noexcept { return superSize[0] * superSize[1] * superSize[2]; }
        /* Static members */
        [[nodiscard]] static EigenSolverType diagonalize(const MatrixType& dynamicMatrix);
        [[nodiscard]] static QPointGrid diagonalize(const MatrixGrid& dynamicMatrixes);
    };

    template<class ScalarType, class PosScalarType>
    PhononSolver<ScalarType, PosScalarType>::PhononSolver(MDCellType unitCell_, Index3D superSize_)
            : unitCell(std::move(unitCell_))
            , superSize(superSize_) {}

    template<class ScalarType, class PosScalarType>
    void PhononSolver<ScalarType, PosScalarType>::swap(PhononSolver& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        unitCell.swap(obj.unitCell);
        superSize.swap(obj.superSize);
    }

    template<class ScalarType, class PosScalarType>
    typename PhononSolver<ScalarType, PosScalarType>::MatrixType PhononSolver<ScalarType, PosScalarType>::interpolatePoint(
            Vector3D qPoint, const MatrixGrid& forceConstants) const {
        const ReciprocalCell repCell = unitCell.reciprocal();
        const Vector3D qVector = repCell.getLattice().transpose() * qPoint;
        const size_t unitCellDOF = getUnitCellDOF();
        FFT3D fft(superSize, {1, 1, 1}, PlanFlag::Estimate);
        auto& rSpace = fft.getRSpace();
        auto& kSpace = fft.getKSpace();
        MatrixType result(unitCellDOF, unitCellDOF);
        for (size_t r = 0; r < unitCellDOF; ++r) {
            for (size_t c = 0; c < unitCellDOF; ++c) {
                kSpace.forIndexInGrid([r, c, &kSpace, &forceConstants](Index3D index) {
                    kSpace(index) = forceConstants(index)(r, c);
                });
                fft.invTransform();
                ComplexType elem = 0;
                rSpace.forIndexInGrid([this, r, c, qVector, &rSpace, &elem](Index3D index) {
                    const auto& lattice = unitCell.getLattice();
                    const Index3D rSpaceDim = rSpace.getDim();
                    ScalarType phase = 0;
                    ScalarType coeff = 1;
                    for (unsigned int i = 0; i < Dim; ++i) {
                        const ssize_t index_i = index[i];
                        const ssize_t dim_i = rSpaceDim[i];
                        const bool isDimOdd = dim_i % 2 != 0;
                        const bool isOnWignerSeitzBoundary = index_i == dim_i;
                        const ScalarType factor = ScalarType(index_i > dim_i / 2 ? index_i - dim_i : (index_i));
                        const ScalarType phase_i = qVector * lattice.row(i).asVector() * factor;
                        if (isDimOdd || !isOnWignerSeitzBoundary)
                            phase += phase_i;
                        else
                            coeff *= cos(phase_i);
                    }
                    const auto factor = ComplexType::fromPhase(phase);
                    elem += rSpace(index) * factor * coeff;
                });
                result(r, c) = elem;
            }
        }
        return result;
    }

    template<class ScalarType, class PosScalarType>
    void PhononSolver<ScalarType, PosScalarType>::toDynamicMatrix(MatrixType& forceConstant) const {
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
    void PhononSolver<ScalarType, PosScalarType>::toDynamicMatrix(MatrixGrid& forceConstants) const {
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
    Vector<ScalarType> PhononSolver<ScalarType, PosScalarType>::makeFreq(const EigenSolverType& eigen) const {
        Vector<ScalarType> result(getUnitCellDOF());
        const Vector<ScalarType> eigenvalues = toRealVector(eigen.getEigenvalues());
        for (size_t i = 0; i < getUnitCellDOF(); ++i)
            result[i] = eigenvalues[i].isNegative() ? -1 : 1;
        result = hadamard(result, sqrt(abs(eigenvalues)));
        result *= ScalarType(1 / (2 * M_PI));
        return result;
    }

    template<class ScalarType, class PosScalarType>
    inline Vector<ScalarType> PhononSolver<ScalarType, PosScalarType>::makeFreq(const QPointGrid& qPoints, Index3D qIndex) const {
        return makeFreq(qPoints(qIndex));
    }

    template<class ScalarType, class PosScalarType>
    DenseMatrix<ScalarType> PhononSolver<ScalarType, PosScalarType>::makeEigenVectors(const EigenSolverType& eigen) const {
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
    inline DenseMatrix<ScalarType> PhononSolver<ScalarType, PosScalarType>::makeEigenVectors(
            const QPointGrid& qPoints, Index3D qIndex) const {
        return makeEigenVectors(qPoints(qIndex));
    }

    template<class ScalarType, class PosScalarType>
    typename PhononSolver<ScalarType, PosScalarType>::EigenSolverType
    PhononSolver<ScalarType, PosScalarType>::diagonalize(const MatrixType& dynamicMatrix) {
        const size_t unitCellDOF = dynamicMatrix.getRow();
        auto eigen = EigenSolverType(unitCellDOF);
        eigen.compute(dynamicMatrix, true);
        eigen.sort();
        return eigen;
    }

    template<class ScalarType, class PosScalarType>
    typename PhononSolver<ScalarType, PosScalarType>::QPointGrid PhononSolver<ScalarType, PosScalarType>::diagonalize(
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
}
