/*
 * Copyright 2020-2024 Weibo He.
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

#include <cstdlib>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

namespace Physica::Core {
    template<Scalar T, int type, size_t maxRow, size_t maxColumn>
    class PLUDecomposition {
    public:
        using MatrixType = DenseMatrix<T, type, maxRow, maxColumn>;
    private:
        MatrixType matrix;
        Array<size_t> biasOrder; //TODO: use permutation matrix instead
    public:
        PLUDecomposition() = default;
        explicit PLUDecomposition(MatrixType m);
        PLUDecomposition(const PLUDecomposition& l) = default;
        PLUDecomposition(PLUDecomposition&& l) noexcept = default;
        ~PLUDecomposition() = default;
        /* Operators */
        PLUDecomposition& operator=(PLUDecomposition obj) noexcept;
        /* Operations */
        void compute(MatrixType m);
        void swap(PLUDecomposition& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const MatrixType& getMatrix() const noexcept { return matrix; }
        [[nodiscard]] const Array<size_t>& getBiasOrder() const noexcept { return biasOrder; }
    private:
        void decompositionColumn(size_t col);
    };
}

#include "PLUDecompositionImpl.h"
