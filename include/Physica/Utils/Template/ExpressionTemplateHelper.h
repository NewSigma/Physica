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
#pragma once

namespace Physica::Utils {
    /**
     * ExpressionType is shared by MatrixExpression and VectorExpression, so it is seperated into a single file.
     */
    enum class ExpressionType {
        Minus,
        Add,
        Sub,
        Mul,
        Div,
        Reciprocal,
        Sqrt,
        Abs,
        Square,
        Ln,
        Sin,
        Cos
    };
}
