/*
 * Copyright 2020-2023 Weibo He.
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

#include "Physica/Core/MultiPrecision/Real.h"
#include "Physica/Core/MultiPrecision/Complex.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"
#include "Physica/Core/Math/NumberTheory/NumberTheory.h"

namespace Physica::Core {
    namespace Internal {
        template<Scalar T>
        T factorial(unsigned int x) {
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

        template<Scalar T>
        T doubleFactorial(size_t x) {
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

    template<Scalar T> T lnGamma(const T& s_);

    template<ScalarOption Option>
    inline Real<Option> gamma(const Real<Option>& s);

    template<ScalarOption Option>
    inline Real<Option> beta(const Real<Option>& s1, const Real<Option>& s2);

    template<Scalar T> T gammaP(const T& a_, const T& x_);

    template<Scalar T> T gammaQ(const T& a_, const T& x_);

    template<Scalar T> T bigamma(const T& x, const T& step);

    template<ScalarOption Option>
    Real<Option> erf(const Real<Option>& x);

    template<Scalar T> T erfc(const T& x_);

    template<ScalarOption Option>
    Real<Option> standardNormalDistribution(const Real<Option>& x);

    template<ScalarOption Option>
    Real<Option> besselJ0(const Real<Option>& x);

    template<ScalarOption Option>
    Real<Option> besselJ1(const Real<Option>& x);

    template<ScalarOption Option>
    Real<Option> besselJn(const Integer& n, const Real<Option>& x);

    template<ScalarOption Option>
    Real<Option> besselJ(const Integer& n, const Real<Option>& x);

    template<ScalarOption Option>
    Real<Option> besselY0(const Real<Option>& x);

    template<ScalarOption Option>
    Real<Option> besselY1(const Real<Option>& x);

    template<ScalarOption Option>
    Real<Option> besselYn(const Integer& n, const Real<Option>& x);

    template<ScalarOption Option>
    void besselJn_Yn_dJn_dYn(
            const Real<Option>& n
            , const Real<Option>& x
            , Real<Option>& __restrict Jn
            , Real<Option>& __restrict Yn
            , Real<Option>& __restrict dJn
            , Real<Option>& __restrict dYn);

    template<ScalarOption Option>
    Real<Option> besselJn(const Real<Option>& n, const Real<Option>& x);

    template<ScalarOption Option>
    Real<Option> besseldJn(const Real<Option>& n, const Real<Option>& x);

    template<ScalarOption Option>
    Real<Option> besselYn(const Real<Option>& n, const Real<Option>& x);

    template<ScalarOption Option>
    Real<Option> besseldYn(const Real<Option>& n, const Real<Option>& x);

    template<ScalarOption Option>
    void sphericalBesselJn_Yn_dJn_dYn(const Real<Option>& n
            , const Real<Option>& x
            , Real<Option>& __restrict jn
            , Real<Option>& __restrict yn
            , Real<Option>& __restrict djn
            , Real<Option>& __restrict dyn);

    template<ScalarOption Option>
    Real<Option> sphericalBesselJn(const Real<Option>& n, const Real<Option>& x);

    template<ScalarOption Option>
    Real<Option> sphericalBesseldJn(const Real<Option>& n, const Real<Option>& x);

    template<ScalarOption Option>
    Real<Option> sphericalBesselYn(const Real<Option>& n, const Real<Option>& x);

    template<ScalarOption Option>
    Real<Option> sphericalBesseldYn(const Real<Option>& n, const Real<Option>& x);

    template<ScalarOption Option>
    Real<Option> legendreP(unsigned int l, const Real<Option>& x);

    template<ScalarOption Option>
    Real<Option> legendreP(unsigned int l, unsigned int m, const Real<Option>& x);

    template<Scalar T>
    Complex<T> sphericalHarmomicY(unsigned int l, int m, const T& theta, const T& phi);
    /**
     * This class generates rotation matrix for spherical hamonic functions
     * 
     * Reference:
     * [1] https://github.com/google/spherical-harmonics.git
     */
    template<Matrix M>
    class HamonicRotator final {
        using ScalarType = M::ScalarType;
    private:
        M initialMat; //Optimize: initialMat may be fixed matrix
        M hamonicRotation; //Current hamonic rotation matrix
    public:
        HamonicRotator(const M& axisRotation);
        HamonicRotator(const HamonicRotator&) = delete;
        HamonicRotator(HamonicRotator&&) = delete;
        ~HamonicRotator() = default;
        /* Operators */
        HamonicRotator& operator=(const HamonicRotator&) = delete;
        HamonicRotator&& operator==(HamonicRotator&&) = delete;
        /* Operations */
        void nextHamonicRotation();
        /* Getters */
        M getCurrentRotation() const { return hamonicRotation; }
    private:
        ScalarType getCenteredElement(size_t row, size_t col);
        bool nearByMargin(double actual, double expected);
        ScalarType P(int i, int a, int b, int l);
        ScalarType U(int m, int n, int l);
        ScalarType V(int m, int n, int l);
        ScalarType W(int m, int n, int l);
    };

    template<Scalar T> T hermiteH(unsigned int n, const T& x);

    template<ScalarOption Option>
    Real<Option> incompBeta(const Real<Option>& a, const Real<Option>& b, const Real<Option>& x);

    template<ScalarOption Option>
    inline Real<Option> studentT(size_t n, const Real<Option>& x);

    template<ScalarOption Option>
    inline Real<Option> distributionF(const Real<Option>& v1, const Real<Option>& v2, const Real<Option>& x);
}

#include "SpetialFunctionsImpl/Bessel.h"
#include "SpetialFunctionsImpl/Gamma.h"
#include "SpetialFunctionsImpl/Legendre.h"
#include "SpetialFunctionsImpl/Hermite.h"
#include "SpetialFunctionsImpl/IncompBeta.h"
