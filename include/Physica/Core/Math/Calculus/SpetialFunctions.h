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
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> lnGamma(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> gamma(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> beta(const Scalar<type, errorTrack>& s1, const Scalar<type, errorTrack>& s2);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> gammaP(const Scalar<type, errorTrack>& a, const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> gammaQ(const Scalar<type, errorTrack>& a, const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> erf(const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> erfc(const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> standardNormalDistribution(const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselJ0(const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselJ1(const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselJn(const Integer& n, const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselJ(const Integer& n, const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselY0(const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselY1(const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselYn(const Integer& n, const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    void besselJn_Yn_dJn_dYn(
            const Scalar<type, errorTrack>& n
            , const Scalar<type, errorTrack>& x
            , Scalar<type, errorTrack>& __restrict Jn
            , Scalar<type, errorTrack>& __restrict Yn
            , Scalar<type, errorTrack>& __restrict dJn
            , Scalar<type, errorTrack>& __restrict dYn);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselJn(const Scalar<type, errorTrack>& n, const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besseldJn(const Scalar<type, errorTrack>& n, const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselYn(const Scalar<type, errorTrack>& n, const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besseldYn(const Scalar<type, errorTrack>& n, const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    void sphericalBesselJn_Yn_dJn_dYn(const Scalar<type, errorTrack>& n
            , const Scalar<type, errorTrack>& x
            , Scalar<type, errorTrack>& __restrict jn
            , Scalar<type, errorTrack>& __restrict yn
            , Scalar<type, errorTrack>& __restrict djn
            , Scalar<type, errorTrack>& __restrict dyn);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> sphericalBesselJn(const Scalar<type, errorTrack>& n, const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> sphericalBesseldJn(const Scalar<type, errorTrack>& n, const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> sphericalBesselYn(const Scalar<type, errorTrack>& n, const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> sphericalBesseldYn(const Scalar<type, errorTrack>& n, const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> legendreP(unsigned int l, const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> legendreP(unsigned int l, unsigned int m, const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    ComplexScalar<type, errorTrack> sphericalHarmomicY(unsigned int l,
                                                int m,
                                                const Scalar<type, errorTrack>& theta,
                                                const Scalar<type, errorTrack>& phi);
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
