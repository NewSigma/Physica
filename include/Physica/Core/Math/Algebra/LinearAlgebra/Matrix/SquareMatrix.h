/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

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
#ifndef PHYSICA_SQUAREMATRIX_H
#define PHYSICA_SQUAREMATRIX_H

#include "Matrix.h"

namespace Physica::Core {
    template<class T, MatrixType type, size_t maxSize>
    class SquareMatrix : public Matrix<T, type, maxSize, maxSize> {
        typedef Matrix<T, type, maxSize, maxSize> Base;
    public:
        enum DeterminateMethod {
            GaussMethod,
            LUMethod
        };
    public:
        SquareMatrix();
        explicit SquareMatrix(size_t length);
        SquareMatrix(const SquareMatrix& m) = default;
        SquareMatrix(SquareMatrix&& m) noexcept;
        ~SquareMatrix() = default;
        /* Operators */
        SquareMatrix& operator=(const SquareMatrix& m) = default;
        SquareMatrix& operator=(SquareMatrix&& m) noexcept;
        /* Getters */
        [[nodiscard]] T determinate(DeterminateMethod method);

        static SquareMatrix getUnitMatrix(size_t length);
    };
}

#include "SquareMatrixImpl.h"

#endif