/*
 * Copyright 2021-2026 Weibo He.
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
using RandomSource = Random<>;

namespace {
    /**
     * All checks are normalized by ||x|| so that they remain meaningful when ||x|| overflows or underflows the naive sum of squares.
     */
    template<Vector V>
    void reflect(const V& x, double prec) noexcept {
        using T = V::ScalarType;
        using Tr = T::RealType;
        const size_t rank = x.getLength();
        V v(rank);
        const Tr norm = x.householder(v);
        const T tau = v[0];
        const T beta = -norm * unit(x[0]);
        v[0] = T(1);

        const V result = x - tau * (v.hermite() * x) * v;
        expect<RandomSource>(scalarNear(result[0] / norm, beta / norm, prec));
        for (size_t i = 1; i < rank; ++i)
            expect<RandomSource>(abs(result[i]) <= Tr(prec) * norm);
    }

    void empty() noexcept {
        using T = float64;
        using VectorType = Vector4D<T>;
        VectorType x{0, 0, 0, 0};
        x.householder();
        for (auto& elem : x)
            expect(elem.isZero());
    }

    void empty_complex() noexcept {
        using T = float64;
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

    void apply() noexcept {
        using T = float64;
        using VectorType = Vector4D<T>;
        const VectorType x{2, 3, 4, 5};
        const size_t rank = x.getLength();
        VectorType v(rank);
        std::ignore = x.householder(v);

        using MatrixType = DenseMatrix<T, MatrixMajor::Col, 4, 4>;
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

    void apply_complex() noexcept {
        using T = float64;
        using Tc = Complex<T>;
        using VectorType = Vector2D<Tc>;
        using MatrixType = DenseMatrix<Tc, MatrixMajor::Col, 2, 2>;

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
    using T = float64;
    // Real
    reflect(Vector4D<T>{2, 3, 4, 5}, 1E-14); // In debug mode, precision can reach 10^-15
    reflect(VectorND<T>::random_uniform<RandomSource>(32), 1E-14);
    reflect(VectorND<Complex<T>>::random_uniform<RandomSource>(32), 1E-14);
    // Extreme magnitude
    reflect(VectorND<T>{1E-300, 1E-200, 1E-150, 1, 1E150, 1E200, 1E300}, 1E-14);
    // Mixed magnitude
    reflect(Vector4D<T>{T(1E200), T(1E-200), T(1), T(0)}, 1E-14);
    reflect(Vector4D<T>{T(1E-200), T(1E200), T(1E-200), T(1)}, 1E-14);
    // Corner cases
    reflect(Vector4D<T>{5, 0, 0, 0}, 1E-14); // Zero tail
    reflect(Vector4D<T>{0, 3, 4, 5}, 1E-14); // Zero leading element
    reflect(Vector4D<T>{-2, 3, -4, 5}, 1E-14); // Mixed signs

    apply();
    apply_complex();
    empty();
    empty_complex();
    return 0;
}
