/*
 * Copyright 2020-2025 Weibo He.
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
    class PLUDecomp {
        using This = PLUDecomp<T, type, maxRow, maxColumn>;
    public:
        using MatrixType = DenseMatrix<T, type, maxRow, maxColumn>;
    private:
        MatrixType matrix;
        Array<size_t> biasOrder; //TODO: use permutation matrix instead
    public:
        PLUDecomp() = default;
        explicit PLUDecomp(MatrixType m);
        PLUDecomp(const This& l) = default;
        PLUDecomp(This&& l) noexcept = default;
        ~PLUDecomp() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void compute(MatrixType m);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const MatrixType& getMatrix() const noexcept { return matrix; }
        [[nodiscard]] const Array<size_t>& getBiasOrder() const noexcept { return biasOrder; }
    private:
        void decomp_col(size_t col);
    };
}

#include "PLUDecompImpl.h"
