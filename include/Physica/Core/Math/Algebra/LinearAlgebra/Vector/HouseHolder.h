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

namespace Physica::Core {
    /**
     * The first element of \param target will be the factor to construct houseHolder matrix.
     * The other parts of \param target will be essential HouseHolder vector.
     * 
     * \return The norm of \param source
     */
    template<class Vector, class OtherVector>
    typename Vector::ScalarType houseHolder(const Vector& source, OtherVector& target) {
        using ScalarType = typename Vector::ScalarType;
        assert(source.getLength() == target.getLength());
        const ScalarType abs_first = abs(source[0]);
        const ScalarType norm = source.norm();
        ScalarType factor = reciprocal(ScalarType(abs_first + norm));
        if (source[0].isNegative())
            factor.toOpposite();

        VectorBlock block_source(source, 1);
        VectorBlock block_target(target, 1);
        block_target = block_source * factor;
        target[0] = ScalarType(1) + abs_first / norm;
        return norm;
    }

    template<class Matrix, class Vector>
    void applyHouseHolder(const Vector& houseHolder, Matrix& mat) {
        Vector copy = houseHolder;
        copy[0] = 1;
        const auto mat1 = (houseHolder[0] * copy).moveToColMatrix();
        mat -= mat1 * (copy.moveToRowMatrix() * mat);
    }

    template<class Matrix, class Vector>
    void applyHouseHolder(Matrix& mat, const Vector& houseHolder) {
        Vector copy = houseHolder;
        copy[0] = 1;
        const auto mat1 = (houseHolder[0] * copy).moveToRowMatrix();
        mat -= (mat * copy.moveToColMatrix()) * mat1;
    }
}
