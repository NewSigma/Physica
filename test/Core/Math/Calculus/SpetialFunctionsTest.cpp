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
#include <iostream>
#include "Physica/Core/Math/Calculus/SpetialFunctions.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

using namespace Physica::Core;

void testLnGamma() {
    constexpr static int count = 2;
    constexpr static double value[count]{13.7, 0.3};
    constexpr static float floatResult[count]{21.77465, 1.095798};
    constexpr static double doubleResult[count]{21.77464517303463, 1.09579799481807552};
    
    for (int i = 0; i < count; ++i) {
        Scalar<Float, false> s(value[i]);
        auto temp = lnGamma(s);
        if ((fabs(float(temp) - floatResult[i]) / fabs(floatResult[i]) >= 1E-6F))
            exit(EXIT_FAILURE);
    }

    for (int i = 0; i < count; ++i) {
        Scalar<Double, false> s(value[i]);
        auto temp = lnGamma(s);
        if ((fabs(double(temp) - doubleResult[i]) / fabs(doubleResult[i]) >= 1E-14))
            exit(EXIT_FAILURE);
    }
}

void testGammaPQ() {
    using T = Scalar<Double, false>;
    constexpr static int count = 2;
    constexpr static double a[count]{13.7, 0.3};
    constexpr static double x[count]{2, 8};
    constexpr static double result[count]{5.309424005280372E-8, 0.99997576072630326};

    for (int i = 0; i < count; ++i) {
        auto temp = gammaP(T(a[i]), T(x[i]));
        if ((fabs(double(temp) - result[i]) >= fabs(result[i]) * 1E-14))
            exit(EXIT_FAILURE);
        temp = gammaQ(T(a[i]), T(x[i]));
        if ((fabs(1 - double(temp) - result[i]) >= fabs(result[i]) * 1E-9))
            exit(EXIT_FAILURE);
    }
}

void testBesselJ() {
    using T = Scalar<Double, false>;
    constexpr static int count = 2;
    constexpr static int n[count]{2, 5};
    constexpr static double x[count]{3, 3};
    constexpr static double result[count]{0.48609126058589107691, 0.043028434877047583925};

    for (int i = 0; i < count; ++i) {
        auto temp = besselJn(n[i], T(x[i]));
        if ((fabs(double(temp) - result[i]) >= fabs(result[i]) * 1E-8))
            exit(EXIT_FAILURE);
    }
}

void testBesselY() {
    using T = Scalar<Double, false>;
    constexpr static int count = 2;
    constexpr static int n[count]{2, 5};
    constexpr static double x[count]{3, 3};
    constexpr static double result[count]{-0.16040039348492372968, -1.9059459538286737322};

    for (int i = 0; i < count; ++i) {
        auto temp = besselYn(n[i], T(x[i]));
        if ((fabs(double(temp) - result[i]) >= fabs(result[i]) * 1E-7))
            exit(EXIT_FAILURE);
    }
}

void testBesselJn_Yn_dJn_dYn() {
    using T = Scalar<Double, false>;
    constexpr static int count = 5;
    constexpr static double n[count]{2, 2, 5, 4, 0.5};
    constexpr static double x[count]{1, 3, 3, 2000, 1};
    constexpr static double result_Jn[count]{0.11490348493190048047, 0.48609126058589107691, 0.043028434877047583925, 0.0070328187752780498324, 0.67139670714180309042};
    constexpr static double result_dJn[count]{0.21024361588113255502, 0.014998118135342407654, 0.060320125796199570454, -0.016398371103788126336, 0.095400514447474534312};
    constexpr static double result_Yn[count]{-1.6506826068162543911, -0.16040039348492372968, -1.9059459538286737322, 0.016396645173086209425, -0.43109886801837607952};
    constexpr static double result_dYn[count]{2.5201523923322200656, 0.43160802044841579822, 2.2598937509893167140, 0.0070287057519738781036, 0.88694614115099113018};

    for (int i = 0; i < count; ++i) {
        T Jn, dJn, Yn, dYn;
        besselJn_Yn_dJn_dYn(T(n[i]), T(x[i]), Jn, Yn, dJn, dYn);
        if ((fabs(double(Jn) - result_Jn[i]) >= fabs(result_Jn[i]) * 1E-9))
            exit(EXIT_FAILURE);
        if ((fabs(double(dJn) - result_dJn[i]) >= fabs(result_dJn[i]) * 1E-10))
            exit(EXIT_FAILURE);
        if ((fabs(double(Yn) - result_Yn[i]) >= fabs(result_Yn[i]) * 1E-10))
            exit(EXIT_FAILURE);
        if ((fabs(double(dYn) - result_dYn[i]) >= fabs(result_dYn[i]) * 1E-9))
            exit(EXIT_FAILURE);
    }
}

void testLegendreP() {
    using T = Scalar<Double, false>;
    constexpr static int count = 2;
    constexpr static unsigned l[count]{5, 4};
    constexpr static unsigned m[count]{2, 3};
    constexpr static double theta[count]{0.37, 0.28};
    constexpr static double answer1[count]{0.30514461613750000, 0.1078912000000};
    constexpr static double answer2[count]{-9.880037322750000, -26.0112384000000};

    for (int i = 0; i < count; ++i) {
        auto result1 = legendreP(l[i], T(theta[i]));
        auto result2 = legendreP(l[i], m[i], T(theta[i]));
        if ((fabs(double(result1) - answer1[i]) >= fabs(answer1[i]) * 1E-15))
            exit(EXIT_FAILURE);
        if ((fabs(double(result2) - answer2[i]) >= fabs(answer2[i]) * 1E-15))
            exit(EXIT_FAILURE);
    }
}

void testSphericalHarmomicY() {
    using T = Scalar<Double, false>;
    constexpr static int count = 2;
    constexpr static unsigned l[count]{5, 4};
    constexpr static unsigned m[count]{2, 3};
    constexpr static double theta[count]{0.37, 0.28};
    constexpr static double phi[count]{0.05, 0.74};
    constexpr static double result_real[count]{0.33052482360605048497, 0.015348915260127907403};
    constexpr static double result_imag[count]{0.03316309979261445896, -0.020223918621792591451};

    for (int i = 0; i < count; ++i) {
        auto result = sphericalHarmomicY(l[i], m[i], T(theta[i]), T(phi[i]));
        if ((fabs(double(result.getReal()) - result_real[i]) >= fabs(result_real[i]) * 1E-5))
            exit(EXIT_FAILURE);
        if ((fabs(double(result.getImag()) - result_imag[i]) >= fabs(result_imag[i]) * 1E-5))
            exit(EXIT_FAILURE);
    }
}
/**
 * Reference:
 * [1] https://github.com/google/spherical-harmonics.git
 */
void testHamonicRotator() {
    constexpr double epsilon = 1E-9;
    using T = Scalar<Double, false>;
    using Matrix = DenseMatrix<T, DenseMatrixOption::Row | DenseMatrixOption::Element>;
    Matrix rotation(3, {0.707106781, -0.707106781, 0, 0.707106781, 0.707106781, 0, 0, 0, 1});
    HamonicRotator rotator(rotation);

    double alpha = M_PI / 4.0;
    /* order 1 */ {
        Matrix answer(3, {cos(alpha), 0, sin(alpha),
                          0, 1, 0,
                          -sin(alpha), 0, cos(alpha)});
        if (!matrixNear(rotator.getCurrentRotation(), answer, epsilon))
            exit(EXIT_FAILURE);
    }
    /* order 2 */ {
        rotator.nextHamonicRotation();
        Matrix answer(5, {cos(2 * alpha), 0, 0, 0, sin(2 * alpha),
                          0, cos(alpha), 0, sin(alpha), 0,
                          0, 0, 1, 0, 0,
                          0, -sin(alpha), 0, cos(alpha), 0,
                          -sin(2 * alpha), 0, 0, 0, cos(2 * alpha)});
        if (!matrixNear(rotator.getCurrentRotation(), answer, epsilon))
            exit(EXIT_FAILURE);
    }
    /* order 3 */ {
        rotator.nextHamonicRotation();
        Matrix answer(7, {cos(3 * alpha), 0, 0, 0, 0, 0, sin(3 * alpha),
                          0, cos(2 * alpha), 0, 0, 0, sin(2 * alpha), 0,
                          0, 0, cos(alpha), 0, sin(alpha), 0, 0,
                          0, 0, 0, 1, 0, 0, 0,
                          0, 0, -sin(alpha), 0, cos(alpha), 0, 0,
                          0, -sin(2 * alpha), 0, 0, 0, cos(2 * alpha), 0,
                          -sin(3 * alpha), 0, 0, 0, 0, 0, cos(3 * alpha)});
        if (!matrixNear(rotator.getCurrentRotation(), answer, epsilon))
            exit(EXIT_FAILURE);
    }
}

int main() {
    testLnGamma();
    testGammaPQ();
    testBesselJ();
    testBesselY();
    testBesselJn_Yn_dJn_dYn();
    testLegendreP();
    testSphericalHarmomicY();
    testHamonicRotator();
    return 0;
}