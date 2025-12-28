/*
 * Copyright 2024-2025 Weibo He.
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

namespace Physica {
    template<Scalar T>
    class PhononPDOS : public PhononDOS<T> {
        using Base = PhononDOS<T>;
        using This = PhononPDOS<T>;
        using typename Base::SolverType;
        using typename Base::MDCellType;
        using typename Base::CoeffVector;
        using typename Base::KSpaceFCGrid;
        using typename Base::EigenValueGrid;
        using Base::Dim;
        using Base::ElementVolume;

        using Base::solver;
        using Base::eigenvalues;
        VectorND<T> direction;
        EigenValueGrid projections;
    public:
        PhononPDOS(MDCellType unitCell, Index3D superSize, const KSpaceFCGrid& forceConstants, Index3D gridDim, VectorND<T> direction_);
        PhononPDOS(const This&) = default;
        PhononPDOS(This&&) noexcept = default;
        ~PhononPDOS() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] T calcPDOS(T freq) const;
        [[nodiscard]] T calcPDOS(T freq, size_t band) const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Base::getUnitCellDOF;
        [[nodiscard]] const auto& getDirection() const noexcept { return direction; }
    };

    template<Scalar T>
    PhononPDOS<T>::PhononPDOS(MDCellType unitCell,
                                       Index3D superSize,
                                       const KSpaceFCGrid& forceConstants,
                                       Index3D gridDim,
                                       VectorND<T> direction_)
            : Base(std::move(unitCell), std::move(superSize), gridDim)
            , direction(std::move(direction_))
            , projections(gridDim) {
        assert(direction.getLength() == getUnitCellDOF() && "[Error]: DOF of direction and unit cell do not match");
        eigenvalues.forND([this, &forceConstants](VectorND<T>& eig, Index3D index) {
            const Index3D shape = eigenvalues.getShape();
            Vector3D<T> qPoint{};
            for (unsigned int i = 0; i < Dim; ++i)
                qPoint[i] = T(index[i]) / T(shape[i]);
            auto fcMatrix = solver.interpolatePoint(qPoint, forceConstants);
            solver.toDynamicMatrix(fcMatrix);
            const auto eigen = SolverType::diagonalize(fcMatrix);
            eig = solver.makeFreq(eigen);

            const size_t unitCellDOF = getUnitCellDOF();
            const auto eigenvectors = solver.makeEigenVectors(eigen);
            VectorND<T> project(unitCellDOF);
            for (size_t i = 0; i < unitCellDOF; ++i)
                project[i] = hadamard(eigenvectors.col(i), direction).squaredNorm();
            projections[index] = project;
        });
    }

    template<Scalar T>
    T PhononPDOS<T>::calcPDOS(T freq) const {
        T result = 0;
        for (size_t band = 0; band < solver.getNumBand(); ++band)
            result += calcPDOS(freq, band);
        return result;
    }

    template<Scalar T>
    T PhononPDOS<T>::calcPDOS(T freq, size_t band) const {
        T result = 0;
        forND(eigenvalues.getShape(), [this, freq, band, &result](Index3D x) {
            const Index3D gridDim = eigenvalues.getShape();
            Index3D y{};
            for (unsigned int i = 0; i < Dim; ++i)
                y[i] = (x[i] + 1) % gridDim[i];

            CoeffVector cornerProj{};
            cornerProj[0] = projections[x][band];
            cornerProj[1] = projections[x[0], x[1], y[2]][band];
            cornerProj[2] = projections[x[0], y[1], x[2]][band];
            cornerProj[3] = projections[x[0], y[1], y[2]][band];
            cornerProj[4] = projections[y[0], x[1], x[2]][band];
            cornerProj[5] = projections[y[0], x[1], y[2]][band];
            cornerProj[6] = projections[y[0], y[1], x[2]][band];
            cornerProj[7] = projections[y][band];
            result += Base::calcElemDOS(freq, band, x) * cornerProj.mean();
        });
        result /= T(eigenvalues.getSize() * ElementVolume);
        return result;
    }

    template<Scalar T>
    void PhononPDOS<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        direction.swap(obj.direction);
        projections.swap(obj.projections);
    }
}
