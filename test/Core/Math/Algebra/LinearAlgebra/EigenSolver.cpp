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
#include "TestHelper.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Transpose.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/EigenSolver.h"

using namespace Physica::Core;

int main() {
    {
        using MatrixType = DenseMatrix<Scalar<Double, false>, DenseMatrixOption::Column | DenseMatrixOption::Vector, 3, 3>;
        const MatrixType mat{{-0.590316, -2.19514, -2.37463},
                            {-1.25006, -0.297493, 1.40349},
                            {0.517063, -0.956614, -0.920775}};
        EigenSolver solver = EigenSolver(mat, true);
        typename EigenSolver<MatrixType>::EigenvalueVector eigenvalues{{-2.13323, 0}, {0.162325, 1.21888}, {0.162325, -1.21888}};
        typename EigenSolver<MatrixType>::EigenvectorMatrix eigenvectors = {{{-0.573583, 0}, {-0.792956, 0}, {-0.205485, 0}},
                                                                            {{-0.225315, -0.221919}, {0.266379, 0.353311}, {0.839165, 0}},
                                                                            {{-0.225315, 0.221919}, {0.266379, -0.353311}, {0.839165, 0}}};
        
    }
    return 0;
}
