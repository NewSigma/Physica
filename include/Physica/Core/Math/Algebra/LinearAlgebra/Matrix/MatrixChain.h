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
#ifndef PHYSICA_MATRIXCHAIN_H
#define PHYSICA_MATRIXCHAIN_H

#include <cstdlib>
#include <qglobal.h>
#include <memory>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"

namespace Physica::Core {
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    class MatrixChain {
        Matrix<T, type, maxRow, maxColumn>** chain;
        size_t** price;
        size_t** point;
        size_t length;
    public:
        explicit MatrixChain(size_t length);
        ~MatrixChain();

        Matrix<T, type, maxRow, maxColumn>*& operator[](size_t i) { Q_ASSERT(i < length); return chain[i]; }
        Matrix<T, type, Dynamic, Dynamic> solve();
    private:
        Matrix<T, type, Dynamic, Dynamic> multiply(size_t from, size_t to);
    };
}

#endif
