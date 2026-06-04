/*
 * Copyright 2020-2026 Weibo He.
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

#include "Physica/Core/Scalar/Real.h"
#include "Physica/Core/Scalar/Complex.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"
#include "Physica/Core/Math/NumberTheory/NumberTheory.h"

namespace Physica {
    namespace Internal {
        template<Scalar T>
        T factorial(unsigned int x) noexcept {
            constexpr static std::array<float64, 16> Cache{1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200, 1307674368000};
            constexpr static int size = Cache.size();
            if (x < size)
                return Cache[x];
            else {
                float64 s = Cache[size - 1];
                float64 i_d = size;
                for (unsigned int i = size; i <= x; i++) {
                    s *= i_d;
                    i_d += float64(1.0);
                }
                return s;
            }
        }

        template<Scalar T>
        T doubleFactorial(size_t x) noexcept {
            constexpr static std::array<float64, 16> Cache{1, 1, 2, 3, 8, 15, 48, 105, 384, 945, 3840, 10395, 46080, 135135, 645120, 2027025};
            constexpr static int size = Cache.size();
            if (x < size)
                return Cache[x];
            else {
                float64 s = Cache[size - 1];
                float64 i_d = size;
                for (size_t i = size; i <= x; i += 2) {
                    s *= i_d;
                    i_d += T(2.0);
                }
                return s;
            }
        }
    }

    template<FloatPrec Prec>
    Real<Prec> lnGamma(const Real<Prec>& x) noexcept;

    template<FloatPrec Prec>
    Real<Prec> gamma(const Real<Prec>& s);

    template<FloatPrec Prec>
    Real<Prec> beta(const Real<Prec>& s1, const Real<Prec>& s2);

    template<Scalar T> T gammaP(const T& a_, const T& x_);

    template<Scalar T> T gammaQ(const T& a_, const T& x_);

    template<FloatPrec Prec>
    Real<Prec> bigamma(const Real<Prec>& x, const Real<Prec>& step) noexcept;

    template<FloatPrec Prec>
    Real<Prec> erf(const Real<Prec>& x);

    template<Scalar T> T erfc(const T& x_);

    template<FloatPrec Prec>
    Real<Prec> standardNormalDistribution(const Real<Prec>& x);

    template<FloatPrec Prec>
    Real<Prec> legendreP(unsigned int l, const Real<Prec>& x);

    template<FloatPrec Prec>
    Real<Prec> legendreP(unsigned int l, unsigned int m, const Real<Prec>& x);

    template<Scalar T>
    Complex<T> sphericalHarmomicY(unsigned int l, int m, const T& theta, const T& phi);
    /**
     * This class generates rotation matrix for spherical harmonic functions
     * 
     * Reference:
     * [1] https://github.com/google/spherical-harmonics.git
     */
    template<Matrix M>
    class HarmonicRotator final {
        using ScalarType = M::ScalarType;
    private:
        M initialMat; //Optimize: initialMat may be fixed matrix
        M harmonicRotation; //Current harmonic rotation matrix
    public:
        HarmonicRotator(const M& axisRotation);
        HarmonicRotator(const HarmonicRotator&) = delete;
        HarmonicRotator(HarmonicRotator&&) = delete;
        ~HarmonicRotator() = default;
        /* Operators */
        HarmonicRotator& operator=(const HarmonicRotator&) = delete;
        HarmonicRotator&& operator==(HarmonicRotator&&) = delete;
        /* Operations */
        void nextHarmonicRotation();
        /* Getters */
        M getCurrentRotation() const { return harmonicRotation; }
    private:
        ScalarType getCenteredElement(size_t row, size_t col);
        bool nearByMargin(double actual, double expected);
        ScalarType P(int i, int a, int b, int l);
        ScalarType U(int m, int n, int l);
        ScalarType V(int m, int n, int l);
        ScalarType W(int m, int n, int l);
    };

    template<Scalar T> T hermiteH(unsigned int n, const T& x);

    template<FloatPrec Prec>
    Real<Prec> incompBeta(const Real<Prec>& a, const Real<Prec>& b, const Real<Prec>& x);

    template<FloatPrec Prec>
    Real<Prec> studentT(size_t n, const Real<Prec>& x);

    template<FloatPrec Prec>
    Real<Prec> distributionF(const Real<Prec>& v1, const Real<Prec>& v2, const Real<Prec>& x);
}

#include "SpecialFunctionsImpl/Bessel.h"
#include "SpecialFunctionsImpl/Gamma.h"
#include "SpecialFunctionsImpl/Legendre.h"
#include "SpecialFunctionsImpl/Hermite.h"
#include "SpecialFunctionsImpl/IncompBeta.h"
