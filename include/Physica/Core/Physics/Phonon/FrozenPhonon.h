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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Grid/RSpaceGrid.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Grid/PeriodIndex3D.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/EigenSolver.h"
#include "Physica/Core/Math/Optimization/OptimizationImpl/QuadraticSearch.h"
#include "Physica/Core/Math/Statistics/NumCharacter.h"
#include "Physica/Core/Math/Transform/FFT.h"
#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"
#include "PhononSolver.h"

namespace Physica::Core {
    /**
     * References:
     * [1] Dario Alfè PHON: A program to calculate phonons using the small displacement method [J]. Computer Physics Communications, 2009, 180(12), 2622-2633 (DOI: 10.1016/j.cpc.2009.03.010)
     */
    template<class ScalarType, class PosScalarType>
    class FrozenPhonon : public PhononSolver<ScalarType, PosScalarType> {
        using This = FrozenPhonon<ScalarType, PosScalarType>;
        using Base = PhononSolver<ScalarType, PosScalarType>;
    public:
        using typename Base::ComplexType;
        using typename Base::Vector3D;
        using typename Base::Index3D;
        using typename Base::FFT3D;
        using typename Base::MDCellType;
        using typename Base::PositionMatrix;
        using typename Base::MatrixType;
        using typename Base::MatrixGrid;
        constexpr static unsigned int Dim = Base::Dim;

        ScalarType displace;
    public:
        FrozenPhonon(MDCellType unitCell, Index3D superSize, ScalarType displace_);
        FrozenPhonon(const FrozenPhonon&) = default;
        FrozenPhonon(FrozenPhonon&&) noexcept = default;
        ~FrozenPhonon() = default;
        /* Operators */
        FrozenPhonon& operator=(FrozenPhonon obj) noexcept;
        /* Operations */
        template<class ForceModel>
        [[nodiscard]] MatrixGrid makeForceConstants(const ForceModel& model) const;
        template<class ForceModel>
        void optimize(const ForceModel& unitCellModel,
                      const ForceModel& superCellModel,
                      ScalarType minFreq,
                      unsigned int numInterpolateStep,
                      unsigned int maxNumStep);
        void swap(FrozenPhonon& obj) noexcept;
        /* Getters */
        using Base::getUnitCell;
        using Base::getNumUnitCellAtom;
        using Base::getUnitCellDOF;
        using Base::getSuperCellDOF;
        using Base::getSuperSize;
        using Base::getForceConstantsGridSize;
        using Base::getNumCell;
    private:
        PositionMatrix makeWignerSeitzRadius() const;
        GridStorage<DenseMatrix<ScalarType>> makeWignerSeitzWeights() const;
        ScalarType calcWignerSeitzWeight(const Vector3D r, const PositionMatrix& wignerSeitzRadius) const;
    };

    template<class ScalarType, class PosScalarType>
    FrozenPhonon<ScalarType, PosScalarType>::FrozenPhonon(MDCellType unitCell, Index3D superSize, ScalarType displace_)
            : Base(std::move(unitCell), superSize), displace(displace_) {}

    template<class ScalarType, class PosScalarType>
    FrozenPhonon<ScalarType, PosScalarType>& FrozenPhonon<ScalarType, PosScalarType>::operator=(FrozenPhonon obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, class PosScalarType>
    template<class ForceModel>
    typename FrozenPhonon<ScalarType, PosScalarType>::MatrixGrid
    FrozenPhonon<ScalarType, PosScalarType>::makeForceConstants(
            const ForceModel& model) const {
        const size_t unitCellDOF = getUnitCellDOF();
        const ScalarType factor = -reciprocal(displace);
        const Index3D superSize = Base::getSuperSize();
        const MDCellType superCell = Base::getUnitCell().template makeSuperCell<ExtendCellOption::CellMajor>(superSize);

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
        return result;
    }

    template<class ScalarType, class PosScalarType>
    template<class ForceModel>
    void FrozenPhonon<ScalarType, PosScalarType>::optimize(
            const ForceModel& unitCellModel,
            const ForceModel& superCellModel,
            ScalarType minFreq,
            unsigned int numInterpolateStep,
            unsigned int maxNumStep) {
        assert(minFreq.isPositive() && "[Error]: Min frequency must be positive");
        for (unsigned int step = 0; step < maxNumStep; ++step) {
            const auto fcMatrixes = makeForceConstants(superCellModel);
            auto fcMatrix = fcMatrixes(0, 0, 0);
            Base::toDynamicMatrix(fcMatrix);
            const auto eigen = Base::diagonalize(fcMatrix);
            const ScalarType freq = Base::makeFreq(eigen)[0];
            const bool isImagFreqSmallEnough = -freq < minFreq;
            if (isImagFreqSmallEnough)
                return;

            const auto eigenvectors = Base::makeEigenVectors(eigen);
            const ScalarType e0 = unitCellModel.potentialEnergy(Base::getUnitCell());
            const ScalarType e1 = unitCellModel.potentialEnergy(Base::shiftAtom(eigenvectors.col(0), displace));
            const ScalarType e2 = unitCellModel.potentialEnergy(Base::shiftAtom(eigenvectors.col(0), -displace));
            QuadraticSearch<ScalarType> optimizer({0, displace, -displace}, {e0, e1, e2});
            for (unsigned int i = 0; i < numInterpolateStep; ++i) {
                optimizer.step([this, &unitCellModel, &eigenvectors](ScalarType dist) {
                    return unitCellModel.potentialEnergy(Base::shiftAtom(eigenvectors.col(0), dist));
                });
            }
            Base::setUnitCell(Base::shiftAtom(eigenvectors.col(0), optimizer.getOptimizedX()));
        }
    }

    template<class ScalarType, class PosScalarType>
    void FrozenPhonon<ScalarType, PosScalarType>::swap(This& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        displace.swap(obj.displace);
    }

    template<class ScalarType, class PosScalarType>
    typename FrozenPhonon<ScalarType, PosScalarType>::PositionMatrix
    FrozenPhonon<ScalarType, PosScalarType>::makeWignerSeitzRadius() const {
        constexpr int OneSideDim = 2; // 2 is enough to generate a Wigner-Seitz cell
        constexpr int TwoSideDim = 2 * OneSideDim + 1;
        constexpr size_t ResultSize = TwoSideDim * TwoSideDim * TwoSideDim - 1;
        const Index3D superSize = Base::getSuperSize();
        PositionMatrix wignerSeitzRadius(ResultSize, 3);
        size_t index = 0;
        for (int i = -OneSideDim; i <= OneSideDim; ++i) {
            for (int j = -OneSideDim; j <= OneSideDim; ++j) {
                for (int k = -OneSideDim; k <= OneSideDim; ++k) {
                    const Vector3D factor{ScalarType(i * superSize[0]), ScalarType(j * superSize[1]), ScalarType(k * superSize[2])};
                    auto row = wignerSeitzRadius.row(index);
                    row = Base::getUnitCell().getLattice().transpose() * factor;
                    const bool isNotGammaPoint = row.squaredNorm() > std::numeric_limits<ScalarType>::epsilon();
                    index += isNotGammaPoint;
                }
            }
        }
        assert(index == ResultSize && "[Error]: Wrong result");
        return wignerSeitzRadius;
    }

    template<class ScalarType, class PosScalarType>
    GridStorage<DenseMatrix<ScalarType>>
    FrozenPhonon<ScalarType, PosScalarType>::makeWignerSeitzWeights() const {
        const auto wignerSeitzRadius = makeWignerSeitzRadius();
        const Index3D superSize = Base::getSuperSize();
        const Index3D gridDim{4 * superSize[0] + 1, 4 * superSize[1] + 1, 4 * superSize[2] + 1};
        const size_t numAtom = getNumUnitCellAtom();
        GridStorage<DenseMatrix<ScalarType>> result(gridDim, numAtom, numAtom);
        GridBase::forIndexInGrid(gridDim, [this, superSize, numAtom, &result, &wignerSeitzRadius](Index3D index) {
            const auto& unitCell = Base::getUnitCell();
            const Vector3D factor{ScalarType(index[0]) - ScalarType(2 * superSize[0]),
                                  ScalarType(index[1]) - ScalarType(2 * superSize[1]),
                                  ScalarType(index[2]) - ScalarType(2 * superSize[2])};
            const Vector3D r0 = unitCell.getLattice().transpose() * factor;
            auto& mat = result(index);
            for (size_t c = 0; c < result.getColumn(); ++c) {
                for (size_t r = 0; r < result.getRow(); ++r) {
                    const Vector3D r1 = r0 + unitCell.getPos().row(r) - unitCell.getPos().row(c);
                    mat(r, c) = calcWignerSeitzWeight(r1, wignerSeitzRadius);
                }
            }
        });
        return result;
    }

    template<class ScalarType, class PosScalarType>
    ScalarType FrozenPhonon<ScalarType, PosScalarType>::calcWignerSeitzWeight(
            const Vector3D r, const PositionMatrix& wignerSeitzRadius) const {
        constexpr double precision = 1E-5;
        size_t count = 1;
        for (size_t i = 0; i < wignerSeitzRadius.getRow(); ++i) {
            const auto radius = wignerSeitzRadius.row(i);
            const ScalarType dot = r * radius.asVector();
            const ScalarType halfSquaredNorm = radius.squaredNorm() * ScalarType(0.5);
            const bool isOnBoundary = scalarNear(dot, halfSquaredNorm, precision);
            if (isOnBoundary) [[unlikely]]
                count += 1;
            else {
                const bool isOutsideWignerSeitzCell = dot > halfSquaredNorm;
                if (isOutsideWignerSeitzCell)
                    return ScalarType(0);
            }
        }
        return reciprocal(ScalarType(count));
    }
}
