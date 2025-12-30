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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Householder.h"
#include "Physica/Core/Scalar/Complex.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using RandomSource = Random<MT19937, 10000>;

namespace {
    template<Vector V>
    static void reflectTest(const V& x, double prec) noexcept {
        using T = V::ScalarType;
        const size_t rank = x.getLength();
        V v(rank);
        const T norm = x.householder(v);
        const T tau = v[0];
        const T beta = -norm * unit(x[0]);
        v[0] = 1;

        V result = x - tau * (v.hermite() * x) * v;
        expect(scalarNear(result[0], beta, prec));

        for (size_t i = 1; i < rank; ++i)
            expect(scalarNear(result[i], T(0), prec));
    }

    static void emptyVectorTest() noexcept {
        using VectorType = Vector4D<T>;
        VectorType x{0, 0, 0, 0};
        x.householder();
        for (auto& elem : x)
            expect(elem.isZero());
    }

    static void emptyComplexVectorTest() noexcept {
        using ComplexType = Complex<T>;
        using ComplexVector = Vector4D<ComplexType>;
        const ComplexVector x = Vector4D<T>{0, 0, 1, 0};
        ComplexVector v(4);
        const T norm = x.householder(v);
        const ComplexType tau = v[0];
        const ComplexType beta = -norm * unit(x[0]);
        v[0] = ComplexType(1);

        const ComplexVector result = x - tau * (v.hermite() * x) * v;
        expect(scalarNear(result[0], beta, 1E-15));

        for (size_t i = 1; i < result.getLength(); ++i)
            expect(scalarNear(result[i], ComplexType(0), 1E-15));
    }

    static void realApplyTest() noexcept {
        using VectorType = Vector4D<T>;
        const VectorType x{2, 3, 4, 5};
        const size_t rank = x.getLength();
        VectorType v(rank);
        std::ignore = x.householder(v);

        using MatrixType = DenseMatrix<T, MatrixOption::Col, 4, 4>;
        const MatrixType m{x, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
        const MatrixType l_answer{{-7.34849, 0, 0, 0}, {-13.0639, 0.203133, -0.729156, -1.66145}, {-20.6846, 0.473976, -1.70137, -3.87671}, {-28.3052, 0.74482, -2.67357, -6.09197}};
        MatrixType l_result = m;

        applyHouseholder(v, l_result);
        expect(matrixNear(l_result, l_answer, 1E-5));

        const MatrixType r_answer{{-16.3299, -18.2351, -20.1402, -22.0454}, {-0.882225, -0.814514, -0.746803, -0.679092}, {1.15703, 0.913982, 0.67093, 0.427878}, {3.19629, 2.64248, 2.08866, 1.53485}};
        MatrixType r_result = m;
        applyHouseholder(r_result, v);
        expect(matrixNear(r_result, r_answer, 1E-5));
    }

    static void complexApplyTest() noexcept {
        using ScalarType = Complex<T>;
        using VectorType = Vector2D<ScalarType>;
        using MatrixType = DenseMatrix<ScalarType, MatrixOption::Col, 2, 2>;

        const VectorType x{{1, 1}, {3, -5}};
        const size_t rank = x.getLength();
        VectorType v(rank);
        const T norm = x.householder(v);

        const MatrixType m{x, {{-2, 7}, {1, 6}}};
        MatrixType householderMat = MatrixType::identity(2);
        applyHouseholder(v, householderMat);
        MatrixType m1 = householderMat * m;
        MatrixType m2 = m;
        applyHouseholder(v, m2);

        expect(scalarNear(m1[0, 0].norm(), norm, 1E-15));
        expect(scalarNear(m1[1, 0].norm(), T(0), 1E-14));

        m1[1, 0] = m2[1, 0] = T(0);
        expect(matrixNear(m1, m2, 1E-14));
        /* Idempotency */ {
            MatrixType m3 = m;
            applyHouseholder(v, m3);
            applyHouseholder(v, m3);
            expect(matrixNear(m, m3, 1E-14));
        }
    }
}

int main() {
    reflectTest(Vector4D<T>{2, 3, 4, 5}, 1E-14); //In debug mode, precision can reach 10^-15
    reflectTest(VectorND<T>::random_uniform<RandomSource>(32), 1E-14);
    reflectTest(VectorND<Complex<T>>::random_uniform<RandomSource>(32), 1E-14);
    realApplyTest();
    complexApplyTest();
    emptyVectorTest();
    emptyComplexVectorTest();
    return 0;
}
