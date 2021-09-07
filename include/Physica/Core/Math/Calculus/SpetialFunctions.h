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

#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Core/MultiPrecision/ComplexScalar.h"
#include "Physica/Core/Math/NumberTheory/NumberTheory.h"

namespace Physica::Core {
    namespace Internal {
        template<class ScalarType>
        ScalarType factorial(unsigned int x) {
            constexpr static int size = 16;
            static const double cache[size] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200, 1307674368000};
            if (x < size)
                return cache[x];
            else {
                double s = cache[size - 1];
                double i_d = size;
                for (unsigned int i = size; i <= x; i++) {
                    s *= i_d;
                    i_d += 1.0;
                }
                return s;
            }
        }

        template<class ScalarType>
        ScalarType doubleFactorial(size_t x) {
            constexpr static size_t size = 16;
            static const double cache[size] = {1, 1, 2, 3, 8, 15, 48, 105, 384, 945, 3840, 10395, 46080, 135135, 645120, 2027025};
            if (x < size)
                return cache[x];
            else {
                double s = cache[size - 1];
                double i_d = size;
                for (size_t i = size; i <= x; i += 2) {
                    s *= i_d;
                    i_d += 2.0;
                }
                return s;
            }
        }
    }

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> lnGamma(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    inline Scalar<option, errorTrack> gamma(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    inline Scalar<option, errorTrack> beta(const Scalar<option, errorTrack>& s1, const Scalar<option, errorTrack>& s2);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> gammaP(const Scalar<option, errorTrack>& a, const Scalar<option, errorTrack>& x);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> gammaQ(const Scalar<option, errorTrack>& a, const Scalar<option, errorTrack>& x);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> erf(const Scalar<option, errorTrack>& x);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> erfc(const Scalar<option, errorTrack>& x);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> standardNormalDistribution(const Scalar<option, errorTrack>& x);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> besselJ0(const Scalar<option, errorTrack>& x);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> besselJ1(const Scalar<option, errorTrack>& x);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> besselJn(const Integer& n, const Scalar<option, errorTrack>& x);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> besselJ(const Integer& n, const Scalar<option, errorTrack>& x);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> besselY0(const Scalar<option, errorTrack>& x);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> besselY1(const Scalar<option, errorTrack>& x);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> besselYn(const Integer& n, const Scalar<option, errorTrack>& x);

    template<ScalarOption option, bool errorTrack>
    void besselJn_Yn_dJn_dYn(
            const Scalar<option, errorTrack>& n
            , const Scalar<option, errorTrack>& x
            , Scalar<option, errorTrack>& __restrict Jn
            , Scalar<option, errorTrack>& __restrict Yn
            , Scalar<option, errorTrack>& __restrict dJn
            , Scalar<option, errorTrack>& __restrict dYn);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> besselJn(const Scalar<option, errorTrack>& n, const Scalar<option, errorTrack>& x);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> besseldJn(const Scalar<option, errorTrack>& n, const Scalar<option, errorTrack>& x);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> besselYn(const Scalar<option, errorTrack>& n, const Scalar<option, errorTrack>& x);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> besseldYn(const Scalar<option, errorTrack>& n, const Scalar<option, errorTrack>& x);

    template<ScalarOption option, bool errorTrack>
    void sphericalBesselJn_Yn_dJn_dYn(const Scalar<option, errorTrack>& n
            , const Scalar<option, errorTrack>& x
            , Scalar<option, errorTrack>& __restrict jn
            , Scalar<option, errorTrack>& __restrict yn
            , Scalar<option, errorTrack>& __restrict djn
            , Scalar<option, errorTrack>& __restrict dyn);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> sphericalBesselJn(const Scalar<option, errorTrack>& n, const Scalar<option, errorTrack>& x);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> sphericalBesseldJn(const Scalar<option, errorTrack>& n, const Scalar<option, errorTrack>& x);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> sphericalBesselYn(const Scalar<option, errorTrack>& n, const Scalar<option, errorTrack>& x);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> sphericalBesseldYn(const Scalar<option, errorTrack>& n, const Scalar<option, errorTrack>& x);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> legendreP(unsigned int l, const Scalar<option, errorTrack>& x);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> legendreP(unsigned int l, unsigned int m, const Scalar<option, errorTrack>& x);

    template<class ScalarType>
    ComplexScalar<ScalarType> sphericalHarmomicY(unsigned int l,
                                                int m,
                                                const ScalarBase<ScalarType>& theta,
                                                const ScalarBase<ScalarType>& phi);
    /**
     * This class generates rotation matrix for spherical hamonic functions
     * 
     * Reference:
     * [1] https://github.com/google/spherical-harmonics.git
     */
    template<class Matrix>
    class HamonicRotator final {
        using ScalarType = typename Matrix::ScalarType;
    private:
        Matrix initialMat; //Optimize: initialMat may be fixed matrix
        Matrix hamonicRotation; //Current hamonic rotation matrix
    public:
        HamonicRotator(const Matrix& axisRotation);
        HamonicRotator(const HamonicRotator&) = delete;
        HamonicRotator(HamonicRotator&&) = delete;
        ~HamonicRotator() = default;
        /* Operators */
        HamonicRotator& operator=(const HamonicRotator&) = delete;
        HamonicRotator&& operator==(HamonicRotator&&) = delete;
        /* Operations */
        void nextHamonicRotation();
        /* Getters */
        Matrix getCurrentRotation() const { return hamonicRotation; }
    private:
        ScalarType getCenteredElement(size_t row, size_t column);
        bool nearByMargin(double actual, double expected);
        ScalarType P(int i, int a, int b, int l);
        ScalarType U(int m, int n, int l);
        ScalarType V(int m, int n, int l);
        ScalarType W(int m, int n, int l);
    };
}

#include "SpetialFunctionsImpl/Bessel.h"
#include "SpetialFunctionsImpl/Gamma.h"
#include "SpetialFunctionsImpl/Legendre.h"
