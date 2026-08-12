/*
 * Copyright 2022-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/LValueMatrix.h"

namespace Physica {
    /**
     * Hold public part of deflatable algorithms. E.g. \class Schur and \class SymmEigenSolver
     */
    class Deflatable {
    protected:
        constexpr static size_t MaxIterationPerCol = 40; // Reference to Eigen

        static size_t activeWindowDownDiag(Matrix auto& mat, size_t upper);
        static size_t activeWindowUpDiag(Matrix auto& mat, size_t upper);
    };
    /**
     * We should process columns whose index is less or equal than \param upper
     *
     * \returns We should process columns whose index is greater or equal to the returned index
     */
    size_t Deflatable::activeWindowDownDiag(Matrix auto& mat, size_t upper) {
        using Trv = std::remove_cvref_t<decltype(mat)>::ScalarType::RealType::ValueType;
        const Trv epsilon = std::numeric_limits<Trv>::epsilon();
        assert(upper < mat.getRow());
        size_t lower = upper;
        size_t lower_1 = upper - 1;
        for (; lower_1 < lower; --lower, --lower_1) { // Make use of overflow
            Trv temp = abs(mat[lower, lower].value()) + abs(mat[lower_1, lower_1].value());
            temp = std::max(temp, epsilon) * epsilon;
            if (abs(mat[lower, lower_1].value()) < temp) {
                mat[lower, lower_1] = 0;
                break;
            }
        }
        return lower;
    }

    size_t Deflatable::activeWindowUpDiag(Matrix auto& mat, size_t upper) {
        using Trv = std::remove_cvref_t<decltype(mat)>::ScalarType::RealType::ValueType;
        assert(upper < mat.getRow());
        size_t lower = upper;
        size_t lower_1 = upper - 1;
        for (; lower_1 < lower; --lower, --lower_1) { // Make use of overflow
            Trv temp = abs(mat[lower, lower].value()) + abs(mat[lower_1, lower_1].value());
            temp = std::max(abs(temp * Trv(std::numeric_limits<Trv>::epsilon())), Trv(std::numeric_limits<Trv>::min()));
            if (abs(mat[lower_1, lower]) < temp) {
                mat[lower_1, lower] = 0;
                break;
            }
        }
        return lower;
    }
}
