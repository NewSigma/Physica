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
#ifndef PHYSICA_LUDECOMPOSITION_H
#define PHYSICA_LUDECOMPOSITION_H

#include <cstdlib>
#include "Matrix/DenseMatrix.h"

namespace Physica::Core {
    template<class T, int type, size_t maxRow, size_t maxColumn>
    class LUDecomposition {
        DenseMatrix<T, type, maxRow, maxColumn> matrix;
        size_t* biasOrder;
    public:
        explicit LUDecomposition(const DenseMatrix<T, type, maxRow, maxColumn>& m);
        explicit LUDecomposition(DenseMatrix<T, type, maxRow, maxColumn>&& m) noexcept;
        LUDecomposition(const LUDecomposition& l);
        LUDecomposition(LUDecomposition&& l) noexcept;
        ~LUDecomposition();
        /* Operators */
        LUDecomposition& operator=(const LUDecomposition& l);
        LUDecomposition& operator=(LUDecomposition&& l) noexcept;
        /* Getters */
        [[nodiscard]] DenseMatrix<T, type, maxRow, maxColumn>&& release() { return std::move(matrix); }
        [[nodiscard]] const DenseMatrix<T, type, maxRow, maxColumn>& getMatrix() { return matrix; }
        [[nodiscard]] const size_t* getOrder() { return biasOrder; }
    private:
        void decompositionColumn(size_t column);
    };
}

#include "LUDecompositionImpl.h"

#endif
