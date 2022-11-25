/*
 * Copyright 2022 WeiBo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/RValueMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrixImpl/HalfDenseMatrixStorage.h"
#include "Physica/Core/Math/Transform/FFT.h"
#include "Physica/Core/Math/Statistics/NumCharacter.h"

namespace Physica::Core {
    class PIPhonon final {
        using ScalarType = Scalar<Double, false>;
        using CorrMatrix = Internal::HalfDenseMatrixStorage<Vector<ScalarType>, Dynamic, Dynamic>;
        constexpr static unsigned int Dim = 3;

        FFT<ScalarType, 3> fft;
        size_t numAtomUnitCell;
        size_t superSizeX;
        size_t superSizeY;
        size_t superSizeZ;
        CorrMatrix force_corr;
        CorrMatrix momentum_corr;
        size_t numSample;
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
        friend std::ostream& operator<<(std::ostream& os, const PIPhonon& phonon);
        friend std::istream& operator>>(std::istream& is, PIPhonon& phonon);
        /* Operations */
        template<class MatrixType>
        void sample(const RValueMatrix<MatrixType>& force, const RValueMatrix<MatrixType>& momentum);
        /* Getters */
        [[nodiscard]] size_t getNumAtomUnitCell() const noexcept { return numAtomUnitCell; }
        [[nodiscard]] size_t getUnitCellDOF() const noexcept { return 3 * getNumAtomUnitCell(); }
        [[nodiscard]] size_t getSuperCellDOF() const noexcept { return getUnitCellDOF() * getNumCell(); }
        [[nodiscard]] size_t getNumCell() const noexcept { return superSizeX * superSizeY * superSizeZ; }
        /* Helpers */
        void swap(PIPhonon& obj) noexcept;
    private:
        void toKSpace();
    };

    template<class MatrixType>
    void PIPhonon::sample(const RValueMatrix<MatrixType>& force, const RValueMatrix<MatrixType>& momentum) {
        assert(force.getRow() == momentum.getRow());
        assert(force.getColumn() == momentum.getColumn());
        assert(force.getRow() == getSuperCellDOF());
        const size_t numCell = getNumCell();
        const size_t unitDof = getUnitCellDOF();
        const size_t superDof = getSuperCellDOF();

        Vector<ScalarType> average_force(superDof);
        Vector<ScalarType> average_momentum(superDof);
        for (size_t i = 0; i < superDof; ++i) {
            average_force[i] = mean(force.row(i));
            average_momentum[i] = mean(momentum.row(i));
        }

        for (size_t r = 0; r < unitDof; ++r) {
            for (size_t c = r; c < unitDof; ++c) {
                for (size_t cell = 0; cell < numCell; ++cell) {
                    const size_t offset_r = r * numCell + cell;
                    const ScalarType force_r = average_force[offset_r];
                    const ScalarType momentum_r = average_momentum[offset_r];

                    const size_t offset_c = c * numCell + cell;
                    const ScalarType force_c = average_force[offset_c];
                    const ScalarType momentum_c = average_momentum[offset_c];

                    const size_t offset_corr = force_corr.accessingIndex(r, c);
                    toNextMean(force_corr[offset_corr][cell], numSample, force_r * force_c);
                    toNextMean(momentum_corr[offset_corr][cell], numSample, momentum_r * momentum_c);
                }
            }
        }
        ++numSample;
    }
}
