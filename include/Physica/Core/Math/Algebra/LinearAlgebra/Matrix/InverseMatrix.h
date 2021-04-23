/*
 * Copyright 2021 WeiBo He.
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
#include <cassert>

namespace Physica::Core {
    template<class MatrixIn>
    class InverseMatrix {
        const MatrixIn& matrixIn;
    public:
        InverseMatrix(const MatrixIn& matrixIn_) : matrixIn(matrixIn_) { assert(matrixIn.getRow() == matrixIn.getColumn()); }
        /* Getters */
        [[nodiscard]] const MatrixIn& getMatrix() const noexcept { return matrixIn; }
        [[nodiscard]] size_t getOrder() const noexcept { return matrixIn.getRow(); }
    };
}