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
#include "Physica/Utils/TestHelper.h"
#include "Physica/Core/Physics/Experiment/DimEstimator.h"

using namespace Physica::Core;
using ScalarType = Scalar<Double, false>;

int main() {
    /**
     * References:
     * [1] Grassberger, P. and Procaccia, I.: Measuring the strangeness of strange attractors, Physica, D9 (1983) 189^208.
     */
    /* logistic */ {
        DenseMatrix<ScalarType, DenseMatrixOption::Column | DenseMatrixOption::Vector, Dynamic, 1> data(30000, 1);
        auto col = data.col(0);
        const ScalarType factor = 3.5699456;
        col[0] = 0.5;
        for (size_t i = 1; i < col.getLength(); ++i) {
            const ScalarType temp = col[i - 1];
            col[i] = factor * temp * (ScalarType(1) - temp);
        }
        const ScalarType dimen = DimEstimator::corrDimen(data, Vector<ScalarType, 8>::linspace(0.00001, 0.0001, 8));
        if (!(ScalarType(0.495) <= dimen && dimen <= ScalarType(0.505)))
            return 1;
    }
    /* Henon map */ {
        DenseMatrix<ScalarType, DenseMatrixOption::Column | DenseMatrixOption::Vector, Dynamic, 2> data(15000, 2);
        auto col1 = data.col(0);
        auto col2 = data.col(1);
        const ScalarType factor1 = 1.4;
        const ScalarType factor2 = 0.3;
        col1[0] = 1;
        col2[0] = 0;
        for (size_t i = 1; i < col1.getLength(); ++i) {
            const ScalarType x = col1[i - 1];
            const ScalarType y = col2[i - 1];
            col1[i] = y + 1 - factor1 * square(x);
            col2[i] = factor2 * x;
        }
        constexpr size_t Length = 8;
        using VectorType = Vector<ScalarType, Length>;
        const VectorType r = exp(VectorType::linspace(-10, -3, Length));
        const ScalarType dimen = DimEstimator::corrDimen(data, r);
        if (!(ScalarType(1.2) <= dimen && dimen <= ScalarType(1.22)))
            return 1;
    }
    return 0;
}
