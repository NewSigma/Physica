/*
 * Copyright 2020-2021 WeiBo He.
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
    template<class T, int type, size_t maxRow, size_t maxColumn>
    class PLUDecomposition {
        DenseMatrix<T, type, maxRow, maxColumn> matrix;
        size_t* biasOrder;
    public:
        explicit PLUDecomposition(const DenseMatrix<T, type, maxRow, maxColumn>& m);
        explicit PLUDecomposition(DenseMatrix<T, type, maxRow, maxColumn>&& m) noexcept;
        PLUDecomposition(const PLUDecomposition& l);
        PLUDecomposition(PLUDecomposition&& l) noexcept;
        ~PLUDecomposition();
        /* Operators */
        PLUDecomposition& operator=(const PLUDecomposition& l);
        PLUDecomposition& operator=(PLUDecomposition&& l) noexcept;
        /* Getters */
        [[nodiscard]] DenseMatrix<T, type, maxRow, maxColumn>&& release() { return std::move(matrix); }
        [[nodiscard]] const DenseMatrix<T, type, maxRow, maxColumn>& getMatrix() const noexcept { return matrix; }
        [[nodiscard]] const size_t* getOrder() { return biasOrder; }
    private:
        void decompositionColumn(size_t column);
    };
}

#include "PLUDecompositionImpl.h"
