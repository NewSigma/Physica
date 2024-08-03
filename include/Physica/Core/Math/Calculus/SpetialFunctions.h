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

    template<class ScalarType>
    ScalarType lnGamma(const ScalarBase<ScalarType>& s_);

    template<ScalarOption Option>
    inline Scalar<Option> gamma(const Scalar<Option>& s);

    template<ScalarOption Option>
    inline Scalar<Option> beta(const Scalar<Option>& s1, const Scalar<Option>& s2);

    template<class ScalarType>
    ScalarType gammaP(const ScalarBase<ScalarType>& a_, const ScalarBase<ScalarType>& x_);

    template<class ScalarType>
    ScalarType gammaQ(const ScalarBase<ScalarType>& a_, const ScalarBase<ScalarType>& x_);

    template<class ScalarType>
    ScalarType bigamma(const ScalarBase<ScalarType>& x, const ScalarType& step);

    template<ScalarOption Option>
    Scalar<Option> erf(const Scalar<Option>& x);

    template<class ScalarType>
    ScalarType erfc(const ScalarBase<ScalarType>& x_);

    template<ScalarOption Option>
    Scalar<Option> standardNormalDistribution(const Scalar<Option>& x);

    template<ScalarOption Option>
    Scalar<Option> besselJ0(const Scalar<Option>& x);

    template<ScalarOption Option>
    Scalar<Option> besselJ1(const Scalar<Option>& x);

    template<ScalarOption Option>
    Scalar<Option> besselJn(const Integer& n, const Scalar<Option>& x);

    template<ScalarOption Option>
    Scalar<Option> besselJ(const Integer& n, const Scalar<Option>& x);

    template<ScalarOption Option>
    Scalar<Option> besselY0(const Scalar<Option>& x);

    template<ScalarOption Option>
    Scalar<Option> besselY1(const Scalar<Option>& x);

    template<ScalarOption Option>
    Scalar<Option> besselYn(const Integer& n, const Scalar<Option>& x);

    template<ScalarOption Option>
    void besselJn_Yn_dJn_dYn(
            const Scalar<Option>& n
            , const Scalar<Option>& x
            , Scalar<Option>& __restrict Jn
            , Scalar<Option>& __restrict Yn
            , Scalar<Option>& __restrict dJn
            , Scalar<Option>& __restrict dYn);

    template<ScalarOption Option>
    Scalar<Option> besselJn(const Scalar<Option>& n, const Scalar<Option>& x);

    template<ScalarOption Option>
    Scalar<Option> besseldJn(const Scalar<Option>& n, const Scalar<Option>& x);

    template<ScalarOption Option>
    Scalar<Option> besselYn(const Scalar<Option>& n, const Scalar<Option>& x);

    template<ScalarOption Option>
    Scalar<Option> besseldYn(const Scalar<Option>& n, const Scalar<Option>& x);

    template<ScalarOption Option>
    void sphericalBesselJn_Yn_dJn_dYn(const Scalar<Option>& n
            , const Scalar<Option>& x
            , Scalar<Option>& __restrict jn
            , Scalar<Option>& __restrict yn
            , Scalar<Option>& __restrict djn
            , Scalar<Option>& __restrict dyn);

    template<ScalarOption Option>
    Scalar<Option> sphericalBesselJn(const Scalar<Option>& n, const Scalar<Option>& x);

    template<ScalarOption Option>
    Scalar<Option> sphericalBesseldJn(const Scalar<Option>& n, const Scalar<Option>& x);

    template<ScalarOption Option>
    Scalar<Option> sphericalBesselYn(const Scalar<Option>& n, const Scalar<Option>& x);

    template<ScalarOption Option>
    Scalar<Option> sphericalBesseldYn(const Scalar<Option>& n, const Scalar<Option>& x);

    template<ScalarOption Option>
    Scalar<Option> legendreP(unsigned int l, const Scalar<Option>& x);

    template<ScalarOption Option>
    Scalar<Option> legendreP(unsigned int l, unsigned int m, const Scalar<Option>& x);

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

    template<class ScalarType>
    ScalarType hermiteH(unsigned int n, const ScalarBase<ScalarType>& x);

    template<ScalarOption Option>
    Scalar<Option> incompBeta(const Scalar<Option>& a, const Scalar<Option>& b, const Scalar<Option>& x);

    template<ScalarOption Option>
    inline Scalar<Option> studentT(size_t n, const Scalar<Option>& x);

    template<ScalarOption Option>
    inline Scalar<Option> distributionF(const Scalar<Option>& v1, const Scalar<Option>& v2, const Scalar<Option>& x);
}

#include "SpetialFunctionsImpl/Bessel.h"
#include "SpetialFunctionsImpl/Gamma.h"
#include "SpetialFunctionsImpl/Legendre.h"
#include "SpetialFunctionsImpl/Hermite.h"
#include "SpetialFunctionsImpl/IncompBeta.h"
