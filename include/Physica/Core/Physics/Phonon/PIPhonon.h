/*
 * Copyright 2022-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseHermiteMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/DenseTensor.h"
#include "Physica/Core/Math/Transform/FFT.h"

namespace Physica {
    class PIPhonon final {
        using ScalarType = float64;
        using ComplexType = Complex<ScalarType>;
        using CorrMatrix = SymmArray<VectorND<ScalarType>, Dynamic>;
        using FFT3D = FFT<ScalarType, 3>;
        constexpr static unsigned int Dim = 3;
        constexpr static double ConsiderAsZeroThrehold = 100 * std::numeric_limits<double>::epsilon();

        FFT3D fft;
        size_t numAtomUnitCell;
        size_t superSizeX;
        size_t superSizeY;
        size_t superSizeZ;

        CorrMatrix force_corr;
        CorrMatrix momentum_corr;
        size_t numSample;

        Array<DenseHermiteMatrix<ComplexType>> kSpaceForceCorr;
        Array<DenseHermiteMatrix<ComplexType>> kSpaceMomentumCorr;

        Array<DenseMatrix<ComplexType>> normalModes;
    public:
        PIPhonon(size_t numAtomUnitCell_,
                 size_t superSizeX_,
                 size_t superSizeY_,
                 size_t superSizeZ_);
        PIPhonon(const PIPhonon&) = default;
        PIPhonon(PIPhonon&&) noexcept = default;
        ~PIPhonon() = default;
        /* Operators */
        PIPhonon& operator=(PIPhonon obj) noexcept;
        /* Operations */
        void sample(const Matrix auto& force, const Matrix auto& momentum);
        void compute();
        void swap(PIPhonon& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumAtomUnitCell() const noexcept { return numAtomUnitCell; }
        [[nodiscard]] size_t getUnitCellDOF() const noexcept { return 3 * getNumAtomUnitCell(); }
        [[nodiscard]] size_t getSuperCellDOF() const noexcept { return getUnitCellDOF() * getNumCell(); }
        [[nodiscard]] size_t getNumCell() const noexcept { return superSizeX * superSizeY * superSizeZ; }
        [[nodiscard]] Index3D getSuperSize() const noexcept { return {superSizeX, superSizeY, superSizeZ}; }
    private:
        void toKSpace();
        void applyTranslationInvariance(DenseHermiteMatrix<ComplexType>& target);
    };

    void PIPhonon::sample(const Matrix auto& force, const Matrix auto& momentum) {
        assert(force.getRow() == momentum.getRow());
        assert(force.getCol() == momentum.getCol());
        assert(force.getRow() == getSuperCellDOF());
        const size_t numCell = getNumCell();
        const size_t unitDof = getUnitCellDOF();
        const size_t superDof = getSuperCellDOF();

        VectorND<ScalarType> average_force(superDof);
        VectorND<ScalarType> average_momentum(superDof);
        for (size_t i = 0; i < superDof; ++i) {
            average_force[i] = force.row(i).mean();
            average_momentum[i] = momentum.row(i).mean();
        }

        VectorND<ScalarType> buffer1(numCell);
        VectorND<ScalarType> buffer2(numCell);
        for (size_t r = 0; r < unitDof; ++r) {
            const size_t offset_r = r * numCell;
            for (size_t c = r; c < unitDof; ++c) {
                const size_t offset_c = c * numCell;
                buffer1 = ScalarType(0);
                buffer2 = ScalarType(0);
                forND(getSuperSize(), [&, this, offset_r, offset_c](Index3D cell1) {
                    forND(getSuperSize(), [&, this, cell1, offset_r, offset_c](Index3D cell2) {
                        const auto range = getSuperSize();
                        Index3D delta;
                        for (int i = 0; i < 3; ++i) {
                            ssize_t temp = static_cast<ssize_t>(cell1[i]) - static_cast<ssize_t>(cell2[i]);
                            if (temp < 0)
                                temp += range[i];
                            else if (temp >= static_cast<ssize_t>(range[i]))
                                temp -= range[i];
                            delta[i] = temp;
                        }

                        const size_t cell1_index1d = cell1[0] * superSizeY * superSizeZ + cell1[1] * superSizeZ + cell1[2];
                        const size_t cell2_index1d = cell2[0] * superSizeY * superSizeZ + cell2[1] * superSizeZ + cell2[2];
                        const size_t offset_buffer = delta[0] * superSizeY * superSizeZ + delta[1] * superSizeZ + delta[2];

                        const size_t offset_r1 = offset_r + cell1_index1d;
                        const ScalarType force_r = average_force[offset_r1];
                        const ScalarType momentum_r = average_momentum[offset_r1];
                        const size_t offset_c1 = offset_c + cell2_index1d;
                        const ScalarType force_c = average_force[offset_c1];
                        const ScalarType momentum_c = average_momentum[offset_c1];
                        buffer1[offset_buffer] += force_r * force_c;
                        buffer2[offset_buffer] += momentum_r * momentum_c;
                    });
                });
                const ScalarType factor = reciprocal(ScalarType(numCell));
                buffer1 *= factor;
                buffer2 *= factor;

                const size_t offset_corr = force_corr.toIndex1D(r, c);
                for (size_t cell = 0; cell < numCell; ++cell) {
                    force_corr.asArray()[offset_corr][cell].toNextMean(numSample, buffer1[cell]);
                    momentum_corr.asArray()[offset_corr][cell].toNextMean(numSample, buffer2[cell]);
                }
            }
        }
        ++numSample;
    }
}
