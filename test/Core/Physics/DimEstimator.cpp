/*
 * Copyright 2022-2025 Weibo He.
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
#include "Physica/Core/Physics/Experiment/DimEstimator.h"

using namespace Physica;
using ScalarType = float64;
/**
 * References:
 * [1] Physica D 9, 189-208 (1983); https://doi.org/10.1016/0167-2789(83)90298-1
 * [2] Phys. Rev. Lett. 45, 1175 (1980); https://doi.org/10.1103/PhysRevLett.45.1175
 */
int main() {
    /* logistic */ {
        DenseMatrix<ScalarType, MatrixOption::Col | MatrixOption::Vector, Dynamic, 1> data(30000, 1);
        auto col = data.col(0);
        const ScalarType factor = 3.5699456;
        col[0] = 0.5;
        for (size_t i = 1; i < col.getLength(); ++i) {
            const ScalarType temp = col[i - 1];
            col[i] = factor * temp * (ScalarType(1) - temp);
        }
        const ScalarType dimen = DimEstimator::corrDimen(data, DenseVector<ScalarType, 8>::linspace(0.00001, 0.0001, 8));
        if (!(ScalarType(0.495) <= dimen && dimen <= ScalarType(0.505)))
            return 1;
    }
    /* Henon map */ {
        // On newer architectures, FMA instructions will affect precision, so we use more data instead of 15000 as [1] does.
        DenseMatrix<ScalarType, MatrixOption::Col | MatrixOption::Vector, Dynamic, 2> data(20000, 2);
        auto col1 = data.col(0);
        auto col2 = data.col(1);
        const ScalarType factor1 = 1.4;
        const ScalarType factor2 = 0.3; // There is a typo in [1]. For supplementary information, refer to [2] (specifically, reference [11] in [1]).
        col1[0] = 1;
        col2[0] = 0;
        for (size_t i = 1; i < col1.getLength(); ++i) {
            const ScalarType x = col1[i - 1];
            const ScalarType y = col2[i - 1];
            col1[i] = y + 1 - factor1 * square(x);
            col2[i] = factor2 * x;
        }
        constexpr size_t Length = 8;
        using VectorType = DenseVector<ScalarType, Length>;
        const VectorType r = exp(VectorType::linspace(-10, -3, Length));
        const ScalarType dimen = DimEstimator::corrDimen(data, r);
        if (!(ScalarType(1.2) <= dimen && dimen <= ScalarType(1.22)))
            return 1;
    }
    return 0;
}
