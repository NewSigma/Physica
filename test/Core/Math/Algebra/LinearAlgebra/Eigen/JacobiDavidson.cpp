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
#include <iostream>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/JacobiDavidson.h"

using namespace Physica::Core;
using ScalarType = Scalar<Double, false>;
using ComplexType = ComplexScalar<ScalarType>;

int main() {
    using Matrix4D = DenseMatrix<ComplexType, MatrixOption::Row | MatrixOption::Vector, 4, 4>;
    const Matrix4D data{
        {{0.7015, -0.8314}, {-1.5771, -2.0026}, {-1.3337, -0.0348}, {0.0229, -0.7145}},
        {{-2.0518, -0.9792}, {0.5080,  0.9642}, {1.1275, -0.7982}, {-0.2620,  1.3514}},
        {{-0.3538, -1.1564}, {0.2820,  0.5201}, {0.3502,  1.0187}, {-1.7502, -0.2248}},
        {{-0.8236, -0.5336}, {0.0335, -0.0200}, {-0.2991, -0.1332}, {-0.2857, -0.5890}},
    };
    const Vector<ScalarType> eigen_answer{-3.5532, -1.7061, 1.1174, 6.6899};
    const DenseHermiteMatrix<ComplexType> hermite = data + data.transpose().conjugate();
    JacobiDavidson<ComplexType> jd(4, 1);
    jd.compute(hermite, Vector<ComplexType>{{1, -2}, {0, -1}, {1, 0}, {-1, 1}});
    std::cout << jd.getEigenvalues().format() << std::endl;
    return 0;
}
