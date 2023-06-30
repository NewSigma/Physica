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
#include "Physica/Core/Physics/Container/KSpaceGrid.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/EigenSolver.h"
#include "Physica/Core/Math/Transform/FFT.h"
#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"

namespace Physica::Core {
    /**
     * References:
     * [1] Dario Alfè PHON: A program to calculate phonons using the small displacement method [J]. Computer Physics Communications, 2009, 180(12), 2622-2633
     */
    template<class ScalarType, class PosScalarType>
    class FinitePhonon {
        using This = FinitePhonon<ScalarType, PosScalarType>;
        using Size3D = Utils::Array<size_t, 3>;
        using MatrixType = DenseMatrix<ComplexScalar<ScalarType>>;
        using EigenSolverType = EigenSolver<MatrixType>;
        using QPointGrid = KSpaceGrid<EigenSolverType>;
    public:
        using MDCellType = MDCell<ScalarType, PosScalarType>;
        using PositionMatrix = typename MDCellType::PositionMatrix;
        constexpr static size_t Dim = Internal::Traits<MDCellType>::Dim;
    private:
        MDCellType unitCell;
        Size3D superSize;
        QPointGrid qPoints;
    public:
        FinitePhonon(MDCellType unitCell_, Size3D superSize_);
        FinitePhonon(const FinitePhonon&) = default;
        FinitePhonon(FinitePhonon&&) noexcept = default;
        ~FinitePhonon() = default;
        /* Operators */
        FinitePhonon& operator=(FinitePhonon obj) noexcept;
        /* Operations */
        template<class ForceModel>
        void diagonalize(const ForceModel& model, ScalarType displace);
        template<class ForceModel>
        KSpaceGrid<MatrixType> makeDynamicMatrix(const ForceModel& model, ScalarType displace);
        void swap(FinitePhonon& obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getUnitCellDOF() const noexcept { return 3 * unitCell.getNumParticle(); }
        [[nodiscard]] size_t getSuperCellDOF() const noexcept { return getUnitCellDOF() * getNumCell(); }
        [[nodiscard]] size_t getNumCell() const noexcept { return superSize[0] * superSize[1] * superSize[2]; }
        [[nodiscard]] const QPointGrid& getQPoints() const noexcept { return qPoints; }
    };

    template<class ScalarType, class PosScalarType>
    FinitePhonon<ScalarType, PosScalarType>::FinitePhonon(MDCellType unitCell_, Size3D superSize_)
            : unitCell(std::move(unitCell_))
            , superSize(superSize_)
            , qPoints(superSize_[0], superSize_[1], superSize_[2]) {}

    template<class ScalarType, class PosScalarType>
    FinitePhonon<ScalarType, PosScalarType>& FinitePhonon<ScalarType, PosScalarType>::operator=(FinitePhonon obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, class PosScalarType>
    template<class ForceModel>
    void FinitePhonon<ScalarType, PosScalarType>::diagonalize(const ForceModel& model, ScalarType displace) {
        auto dynamicMatrixGrid = makeDynamicMatrix(model, displace);
        for (size_t i = 0; i < dynamicMatrixGrid.flatten().getLength(); ++i)
            qPoints.flatten()[i].compute(dynamicMatrixGrid.flatten()[i], true);
    }

    template<class ScalarType, class PosScalarType>
    template<class ForceModel>
    KSpaceGrid<typename FinitePhonon<ScalarType, PosScalarType>::MatrixType>
    FinitePhonon<ScalarType, PosScalarType>::makeDynamicMatrix(const ForceModel& model, ScalarType displace) {
        const ScalarType factor = -reciprocal(displace);
        MDCellType superCell = unitCell;
        superCell.template unitToSuper<ExtendCellOption::CellMajor>(superSize[0], superSize[1], superSize[2]);
        PositionMatrix pos = superCell.getPos();
        KSpaceGrid<MatrixType> dynamicMatrixGrid(superSize[0], superSize[1], superSize[2]);
        FFT<ScalarType, 1> fft(getNumCell(), 1);
        for (size_t row = 0; row < getUnitCellDOF(); ++row) {
            ScalarType& toDisplace = pos(row / Dim, row % Dim);
            const ScalarType copy = toDisplace;
            toDisplace += displace;
            const Vector<ScalarType> forceConst =
                    model.template force<SequentialExecutor>(MDCellType(superCell.getLattice(), pos, superCell.getMassVec())) * factor;
            toDisplace = copy;

            for (size_t col = row; col < getUnitCellDOF(); ++col) {
                const size_t shift = getUnitCellDOF();
                for (size_t cell = 0; cell < getNumCell(); ++cell)
                    fft.getRSpace()[cell] = forceConst[col + cell * shift];
                fft.transform();
                auto& matrixArray = dynamicMatrixGrid.flatten();
                for (size_t i = 0; i < matrixArray.getLength(); ++i)
                    matrixArray[i](row, col) = fft.getKSpace()[i];
            }
        }
        return dynamicMatrixGrid;
    }

    template<class ScalarType, class PosScalarType>
    void FinitePhonon<ScalarType, PosScalarType>::swap(This& obj) noexcept {
        unitCell.swap(obj.unitCell);
        superSize.swap(obj.superSize);
        qPoints.swap(obj.qPoints);
    }
}
