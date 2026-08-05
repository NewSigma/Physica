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
#include "Physica/Core/Math/Calculus/SpecialFunctions.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;

namespace {
    void testLnGamma() noexcept {
        constexpr static Array<double, 2> value{13.7, 0.3};
        constexpr static Array<float, 2> floatResult{21.77465, 1.095798};
        constexpr static Array<double, 2> doubleResult{21.77464517303463, 1.09579799481807552};
        for (size_t i = 0; i < value.size(); ++i) {
            using T = float32;
            T s(value[i]);
            auto temp = lnGamma(s);
            expect(scalarNear(temp, T(floatResult[i]), 1E-6));
        }

        for (size_t i = 0; i < value.size(); ++i) {
            using T = float64;
            T s(value[i]);
            auto temp = lnGamma(s);
            expect(scalarNear(temp, T(doubleResult[i]), 1E-15));
        }
    }

    void testGammaPQ() noexcept {
        using T = float64;
        constexpr static Array<double, 2> a{13.7, 0.3};
        constexpr static Array<double, 2> x{2, 8};
        constexpr static Array<double, 2> result{5.309424005280372E-8, 0.99997576072630326};
        for (size_t i = 0; i < a.size(); ++i) {
            auto temp = gammaP(T(a[i]), T(x[i]));
            expect(scalarNear(temp, T(result[i]), 1E-14));

            temp = gammaQ(T(a[i]), T(x[i]));
            expect(scalarNear(temp, T(1 - result[i]), 1E-11));
        }
        // Test special values
        expect<RandomSource>(gammaP(T::random_uniform<RandomSource>(), T(0)).isZero());
        expect<RandomSource>(gammaQ(T::random_uniform<RandomSource>(), T(0)) == T(1));
    }

    void testBiGamma() noexcept {
        using T = float64;
        constexpr static Array<float64, 2> x{13.7, 0.3};
        constexpr static Array<float64, 2> step{1E-1, 1E-1};
        constexpr static Array<float64, 2> result{2.5804557238996526, -3.5025242222001330};
        for (size_t i = 0; i < x.size(); ++i) {
            auto temp = bigamma(T(x[i]), T(step[i]));
            expect(scalarNear(temp, T(result[i]), 1E-13));
        }
    }

    void testBesselJ() noexcept {
        using T = float64;
        constexpr static Array<int, 2> n{2, 5};
        constexpr static Array<double, 2> x{3, 3};
        constexpr static Array<double, 2> result{0.48609126058589107691, 0.043028434877047583925};
        for (size_t i = 0; i < n.size(); ++i) {
            auto temp = besselJn(n[i], T(x[i]));
            expect(scalarNear(temp, T(result[i]), 1E-8));
        }
    }

    void testBesselY() noexcept {
        using T = float64;
        constexpr static Array<int, 2> n{2, 5};
        constexpr static Array<double, 2> x{3, 3};
        constexpr static Array<double, 2> result{-0.16040039348492372968, -1.9059459538286737322};
        for (size_t i = 0; i < n.size(); ++i) {
            auto temp = besselYn(n[i], T(x[i]));
            expect(scalarNear(temp, T(result[i]), 1E-7));
        }
    }

    void testBesselJn_Yn_dJn_dYn() noexcept {
        using T = float64;
        constexpr static Array<double, 5> n{2, 2, 5, 4, 0.5};
        constexpr static Array<double, 5> x{1, 3, 3, 2000, 1};
        constexpr static Array<double, 5> result_Jn{0.11490348493190048047, 0.48609126058589107691, 0.043028434877047583925, 0.0070328187752780498324, 0.67139670714180309042};
        constexpr static Array<double, 5> result_dJn{0.21024361588113255502, 0.014998118135342407654, 0.060320125796199570454, -0.016398371103788126336, 0.095400514447474534312};
        constexpr static Array<double, 5> result_Yn{-1.6506826068162543911, -0.16040039348492372968, -1.9059459538286737322, 0.016396645173086209425, -0.43109886801837607952};
        constexpr static Array<double, 5> result_dYn{2.5201523923322200656, 0.43160802044841579822, 2.2598937509893167140, 0.0070287057519738781036, 0.88694614115099113018};
        for (size_t i = 0; i < n.size(); ++i) {
            T Jn, dJn, Yn, dYn;
            besselJn_Yn_dJn_dYn(T(n[i]), T(x[i]), Jn, Yn, dJn, dYn);
            expect(scalarNear(Jn, T(result_Jn[i]), 1E-9));
            expect(scalarNear(dJn, T(result_dJn[i]), 1E-10));
            expect(scalarNear(Yn, T(result_Yn[i]), 1E-10));
            expect(scalarNear(dYn, T(result_dYn[i]), 1E-9));
        }
    }

    void testLegendreP() noexcept {
        using T = float64;
        constexpr static Array<unsigned, 2> l{5, 4};
        constexpr static Array<unsigned, 2> m{2, 3};
        constexpr static Array<double, 2> theta{0.37, 0.28};
        constexpr static Array<double, 2> answer1{0.30514461613750000, 0.1078912000000};
        constexpr static Array<double, 2> answer2{-9.880037322750000, -26.0112384000000};
        for (size_t i = 0; i < l.size(); ++i) {
            auto result1 = legendreP(l[i], T(theta[i]));
            auto result2 = legendreP(l[i], m[i], T(theta[i]));
            expect(scalarNear(result1, T(answer1[i]), 1E-15));
            expect(scalarNear(result2, T(answer2[i]), 1E-15));
        }
    }

    void testSphericalHarmomicY() noexcept {
        using T = float64;
        constexpr static Array<unsigned, 2> l{5, 4};
        constexpr static Array<int, 2> m{2, 3};
        constexpr static Array<double, 2> theta{0.37, 0.28};
        constexpr static Array<double, 2> phi{0.05, 0.74};
        constexpr static Array<double, 2> result_real{0.33052482360605048497, 0.015348915260127907403};
        constexpr static Array<double, 2> result_imag{0.03316309979261445896, -0.020223918621792591451};
        for (size_t i = 0; i < l.size(); ++i) {
            auto result = sphericalHarmomicY(l[i], m[i], T(theta[i]), T(phi[i]));
            expect(scalarNear(result, Complex<T>(T(result_real[i]), T(result_imag[i])), 1E-15));
        }
    }
    /**
    * Reference:
    * [1] https://github.com/google/spherical-harmonics.git
    */
    void testHarmonicRotator() noexcept {
        constexpr double epsilon = 1E-9;
        using T = float64;
        using MatrixType = DenseMatrix<T, MatrixMajor::Row>;
        using Matrix3D = DenseMatrix<T, MatrixMajor::Row, 3, 3>;
        const Matrix3D rotation({0.707106781, -0.707106781, 0, 0.707106781, 0.707106781, 0, 0, 0, 1});
        HarmonicRotator<MatrixType> rotator(rotation);

        constexpr double alpha = M_PI / 4.0;
        /* order 1 */ {
            const Matrix3D answer({cos(alpha), 0, sin(alpha),
                                0, 1, 0,
                                -sin(alpha), 0, cos(alpha)});
            expect(matrixNear(rotator.getCurrentRotation(), answer, epsilon));
        }
        /* order 2 */ {
            using Matrix5D = DenseMatrix<T, MatrixMajor::Row, 5, 5>;
            rotator.nextHarmonicRotation();
            const Matrix5D answer({cos(2 * alpha), 0, 0, 0, sin(2 * alpha),
                                0, cos(alpha), 0, sin(alpha), 0,
                                0, 0, 1, 0, 0,
                                0, -sin(alpha), 0, cos(alpha), 0,
                                -sin(2 * alpha), 0, 0, 0, cos(2 * alpha)});
            expect(matrixNear(rotator.getCurrentRotation(), answer, epsilon));
        }
        /* order 3 */ {
            using Matrix7D = DenseMatrix<T, MatrixMajor::Row, 7, 7>;
            rotator.nextHarmonicRotation();
            const Matrix7D answer({cos(3 * alpha), 0, 0, 0, 0, 0, sin(3 * alpha),
                                0, cos(2 * alpha), 0, 0, 0, sin(2 * alpha), 0,
                                0, 0, cos(alpha), 0, sin(alpha), 0, 0,
                                0, 0, 0, 1, 0, 0, 0,
                                0, 0, -sin(alpha), 0, cos(alpha), 0, 0,
                                0, -sin(2 * alpha), 0, 0, 0, cos(2 * alpha), 0,
                                -sin(3 * alpha), 0, 0, 0, 0, 0, cos(3 * alpha)});
            expect(matrixNear(rotator.getCurrentRotation(), answer, epsilon));
        }
    }

    void testHermiteH() noexcept {
        using T = float64;
        constexpr static Array<unsigned int, 3> n{2, 5, 3};
        constexpr static Array<double, 3> x{3, 3, 18};
        constexpr static Array<double, 3> result{34, 3816, 46440};
        for (size_t i = 0; i < n.size(); ++i) {
            auto temp = hermiteH(n[i], T(x[i]));
            expect(scalarNear(temp, T(result[i]), 1E-14));
        }
    }

    void testIncompBeta() noexcept {
        using T = float64;
        constexpr static Array<unsigned int, 3> n{5, 10, 12};
        constexpr static Array<double, 3> x{1.2, 0.2, 2.1};
        constexpr static Array<double, 3> result{0.7161089432938979, 0.1545108086196544, 0.94245506126504887};
        for (size_t i = 0; i < n.size(); ++i) {
            auto temp = studentT(n[i], T(x[i]));
            expect(scalarNear(temp, T(result[i]), 1E-13));
        }
    }
}

int main() {
    testLnGamma();
    testGammaPQ();
    testBiGamma();
    testBesselJ();
    testBesselY();
    testBesselJn_Yn_dJn_dYn();
    testLegendreP();
    testSphericalHarmomicY();
    testHarmonicRotator();
    testHermiteH();
    testIncompBeta();
    return 0;
}
