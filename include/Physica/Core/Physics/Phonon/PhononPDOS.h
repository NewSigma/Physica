/*
 * Copyright 2024 Weibo He.
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

#include "PhononDOS.h"

namespace Physica::Core {
    template<class ScalarType>
    class PhononPDOS : public PhononDOS<ScalarType> {
        using Base = PhononDOS<ScalarType>;
        using typename Base::SolverType;
        using typename Base::MDCellType;
        using typename Base::Vector3D;
        using typename Base::CoeffVector;
        using typename Base::Index3D;
        using typename Base::KSpaceFCGrid;
        using typename Base::EigenValueGrid;
        using Base::Dim;
        using Base::ElementVolume;
        using VectorType = VectorND<ScalarType>;

        using Base::solver;
        using Base::eigenvalues;
        VectorType direction;
        EigenValueGrid projections;
    public:
        PhononPDOS(MDCellType unitCell, Index3D superSize, const KSpaceFCGrid& forceConstants, Index3D gridDim, VectorType direction_);
        PhononPDOS(const PhononPDOS&) = default;
        PhononPDOS(PhononPDOS&&) noexcept = default;
        ~PhononPDOS() = default;
        /* Operators */
        PhononPDOS& operator=(PhononPDOS obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] ScalarType calcPDOS(ScalarType freq) const;
        [[nodiscard]] ScalarType calcPDOS(ScalarType freq, size_t band) const;
        void swap(PhononPDOS& __restrict obj) noexcept;
        /* Getters */
        using Base::getUnitCellDOF;
        [[nodiscard]] const VectorType& getDirection() const noexcept { return direction; }
    };

    template<class ScalarType>
    PhononPDOS<ScalarType>::PhononPDOS(MDCellType unitCell,
                                       Index3D superSize,
                                       const KSpaceFCGrid& forceConstants,
                                       Index3D gridDim,
                                       VectorType direction_)
            : Base(std::move(unitCell), std::move(superSize), gridDim)
            , direction(std::move(direction_))
            , projections(gridDim) {
        assert(direction.getLength() == getUnitCellDOF() && "[Error]: DOF of direction and unit cell do not match");
        eigenvalues.forIndexInGrid([this, &forceConstants](Index3D index) {
            const Index3D gridDim = eigenvalues.getDim();
            Vector3D qPoint{};
            for (unsigned int i = 0; i < Dim; ++i)
                qPoint[i] = ScalarType(index[i]) / ScalarType(gridDim[i]);
            auto fcMatrix = solver.interpolatePoint(qPoint, forceConstants);
            solver.toDynamicMatrix(fcMatrix);
            const auto eigen = SolverType::diagonalize(fcMatrix);
            eigenvalues(index) = solver.makeFreq(eigen);

            const size_t unitCellDOF = getUnitCellDOF();
            const auto eigenvectors = solver.makeEigenVectors(eigen);
            VectorType project(unitCellDOF);
            for (size_t i = 0; i < unitCellDOF; ++i)
                project[i] = hadamard(eigenvectors.col(i), direction).squaredNorm();
            projections(index) = project;
        });
    }

    template<class ScalarType>
    ScalarType PhononPDOS<ScalarType>::calcPDOS(ScalarType freq) const {
        ScalarType result = 0;
        for (size_t band = 0; band < solver.getNumBand(); ++band)
            result += calcPDOS(freq, band);
        return result;
    }

    template<class ScalarType>
    ScalarType PhononPDOS<ScalarType>::calcPDOS(ScalarType freq, size_t band) const {
        ScalarType result = 0;
        eigenvalues.forIndexInGrid([this, freq, band, &result](Index3D index) {
            const Index3D gridDim = eigenvalues.getDim();
            Index3D index1{};
            for (unsigned int i = 0; i < Dim; ++i)
                index1[i] = (index[i] + 1) % gridDim[i];

            CoeffVector cornerProj{};
            cornerProj[0] = projections(index)[band];
            cornerProj[1] = projections(index[0], index[1], index1[2])[band];
            cornerProj[2] = projections(index[0], index1[1], index[2])[band];
            cornerProj[3] = projections(index[0], index1[1], index1[2])[band];
            cornerProj[4] = projections(index1[0], index[1], index[2])[band];
            cornerProj[5] = projections(index1[0], index[1], index1[2])[band];
            cornerProj[6] = projections(index1[0], index1[1], index[2])[band];
            cornerProj[7] = projections(index1)[band];
            result += Base::calcElemDOS(freq, band, index) * mean(cornerProj);
        });
        result /= ScalarType(eigenvalues.getSize() * ElementVolume);
        return result;
    }

    template<class ScalarType>
    void PhononPDOS<ScalarType>::swap(PhononPDOS& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        direction.swap(obj.direction);
        projections.swap(obj.projections);
    }
}
