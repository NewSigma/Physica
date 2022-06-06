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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/LValueMatrix.h"

namespace Physica::Core {
    /**
     * Hold public part of decouplable algorithms. E.g. \class Schur and \class SymmEigenSolver
     */
    class Decouplable {
    protected:
        constexpr static size_t maxItePerCol = 40; //Reference to Eigen

        template<class MatrixType>
        static size_t activeWindowDownDiag(LValueMatrix<MatrixType>& __restrict mat, size_t upper);
        template<class MatrixType>
        static size_t activeWindowUpDiag(LValueMatrix<MatrixType>& __restrict mat, size_t upper);
    };

    /**
     * We should process columns whose index is less or equal than \param upper
     * 
     * \returns We should process columns whose index is greater or equal to the returned index
     */
    template<class MatrixType>
    size_t Decouplable::activeWindowDownDiag(LValueMatrix<MatrixType>& __restrict mat, size_t upper) {
        using ScalarType = typename MatrixType::ScalarType;
        using RealType = typename ScalarType::RealType;

        assert(upper < mat.getRow());
        size_t lower = upper;
        size_t lower_1 = upper - 1;
        for (; lower_1 < lower; --lower, --lower_1) { //Make use of overflow
            RealType temp = abs(mat(lower, lower)) + abs(mat(lower_1, lower_1));
            temp = std::max(abs(temp * RealType(std::numeric_limits<ScalarType>::epsilon())), RealType(std::numeric_limits<ScalarType>::min()));
            if (abs(mat(lower, lower_1)) < temp) {
                mat(lower, lower_1) = ScalarType::Zero();
                break;
            }
        }
        return lower;
    }

    template<class MatrixType>
    size_t Decouplable::activeWindowUpDiag(LValueMatrix<MatrixType>& __restrict mat, size_t upper) {
        using ScalarType = typename MatrixType::ScalarType;
        using RealType = typename ScalarType::RealType;

        assert(upper < mat.getRow());
        size_t lower = upper;
        size_t lower_1 = upper - 1;
        for (; lower_1 < lower; --lower, --lower_1) { //Make use of overflow
            RealType temp = abs(mat(lower, lower)) + abs(mat(lower_1, lower_1));
            temp = std::max(abs(temp * RealType(std::numeric_limits<ScalarType>::epsilon())), RealType(std::numeric_limits<ScalarType>::min()));
            if (abs(mat(lower_1, lower)) < temp) {
                mat(lower_1, lower) = ScalarType::Zero();
                break;
            }
        }
        return lower;
    }
}
