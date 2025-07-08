/*
 * Copyright 2021-2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/UnitMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/EigenSolver.h"
#include "Physica/Core/Physics/ManyBody/Hamilton/TransIsingMatrix.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/SpinRepr.h"

using namespace Physica;

template<Vector V>
bool vectorNearZero(const LValueVector<V>& v, double precision) {
    using ScalarType = V::ScalarType;
    for (size_t i = 0; i < v.getLength(); ++i)
        if (!scalarNear(v[i], ScalarType(0), precision))
            return false;
    return true;
}

template<Matrix M, Backend B>
bool eigenTestImpl(const M& mat, double precision) {
    using ScalarType = M::ScalarType;
    using ComplexVector = EigenSolver<ScalarType>::EigenvalueVector;
    auto solver = EigenSolver<ScalarType>(mat.getRow());
    if constexpr (B == Backend::Base)
        solver.compute_base(mat, true);
    else
        solver.compute_mkl(mat, true);

    const size_t order = mat.getRow();
    auto eigenvectors = solver.getEigenvectors();
    for (size_t i = 0; i < order; ++i) {
        const ComplexVector result = (mat - solver.getEigenvalues()[i] * UnitMatrix<ScalarType>(order)) * eigenvectors.col(i);
        if (!vectorNearZero(result, precision))
            return false;
    }

    solver.sort();
    eigenvectors = solver.getEigenvectors();
    for (size_t i = 0; i < order; ++i) {
        if (i > 1 && solver.getEigenvalues()[i - 1].real() > solver.getEigenvalues()[i].real())
            return false;
        const ComplexVector result = (mat - solver.getEigenvalues()[i] * UnitMatrix<ScalarType>(order)) * eigenvectors.col(i);
        if (!vectorNearZero(result, precision))
            return false;
    }
    return true;
}

template<Matrix M>
bool eigenTest(const M& mat, double precision) {
    using ScalarType = M::ScalarType;
    if (!eigenTestImpl<M, Backend::Base>(mat, precision))
        return false;
    if constexpr (HasMKL() && (ScalarType::Prec == Float32 || ScalarType::Prec == Float64)) {
        if (!eigenTestImpl<M, Backend::MKL>(mat, precision))
            return false;
    }
    return true;
}

template<Matrix M, bool IsHermite>
bool reconstructTest(const M& mat, double precision) {
    using ScalarType = M::ScalarType;
    auto solver = EigenSolver<ScalarType>(mat, true);
    const M result = IsHermite ? solver.reconstruct_hermite() : solver.reconstruct();
    if (!matrixNear(result, mat, precision))
        return false;
    return true;
}

int main() {
    using RealType = float64;
    using ComplexType = Complex<RealType>;
    {
        using T = DenseMatrix<RealType, MatrixOption::Col | MatrixOption::Vector, 3, 3>;
        const T mat1{{-0.590316, -2.19514, -2.37463},
                     {-1.25006, -0.297493, 1.40349},
                     {0.517063, -0.956614, -0.920775}};
        if (!eigenTest(mat1, 1E-14))
            return 1;
        if (!reconstructTest<T, false>(mat1, 1E-14))
            return 1;

        const T mat2{{0.354784, 0.604092, 0.557408},
                     {0.221811, 0.484944, 0.181502},
                     {0.491429, 0.0161911, 0.772586}};
        if (!eigenTest(mat2, 1E-15))
            return 1;
        
        const T mat3 = mat2 + mat2.transpose();
        if (!reconstructTest<T, true>(mat3, 1E-14))
            return 1;
    }
    {
        using T = DenseMatrix<RealType, MatrixOption::Col | MatrixOption::Vector, 5, 5>;
        const T mat1{{0.200743, 0.151314, 0.152894, 0.934051, 0.487404},
                     {0.819659, 0.434558, 0.829935, 0.837801, 0.699088},
                     {0.432202, 0.744724, 0.823444, 0.149924, 0.72579},
                     {0.315664, 0.789745, 0.892616, 0.446226, 0.0861119},
                     {0.569445, 0.734035, 0.816045, 0.538242, 0.302492}};
        if (!eigenTest(mat1, 1E-14))
            return 1;

        const T mat2 {{0.502016, 0.812053, 0.0716913, 0.180531, 0.0161382},
                      {0.54811, 0.428973, 0.997804, 0.601252, 0.0704318},
                      {0.885699, 0.0248232, 0.397405, 0.0504956, 0.423087},
                      {0.645844, 0.650145, 0.264407, 0.264864, 0.897732},
                      {0.244945, 0.716587, 0.714982, 0.228202, 0.221696}};
        if (!eigenTest(mat2, 1E-14))
            return 1;
    }
    /* Test degeneracy */ {
        using T = DenseMatrix<RealType, MatrixOption::Col | MatrixOption::Vector, 6, 6>;
        const T mat1{{ 0.1343184046,             0,             0, -0.1343184056,             0,             0},
                     {            0,  0.1341424528,             0,             0, -0.1341424541,             0},
                     {            0,             0,  0.1342191829,             0,             0, -0.1342191848},
                     {-0.1343184056,             0,             0,  0.1343184065,             0,             0},
                     {            0, -0.1341424541,             0,             0,  0.1341424554,             0},
                     {            0,             0, -0.1342191848,             0,             0,  0.1342191868}};
        if (!eigenTest(mat1, 1E-15))
            return 1;

        using ComplexMatrix = DenseMatrix<ComplexType, MatrixOption::Col | MatrixOption::Vector, 6, 6>;
        const ComplexMatrix mat2 = mat1;
        if (!eigenTest(mat2, 1E-15))
            return 1;
    }
    /* Test degeneracy */ {
        DenseMatrix<float64> m = TransIsingMatrix(TransIsing<float64, 1>({{8}, 1}, 1, 0.01), SpinRepr<1, 8>(8));
        if (!eigenTest(m, 1E-13))
            return 1;
    }
    return 0;
}
