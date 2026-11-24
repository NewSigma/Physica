/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Scalar/Real.h"
#include "Physica/Core/Scalar/Complex.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"
#include "Physica/Core/Math/NumberTheory/NumberTheory.h"

namespace Physica {
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
        HarmonicRotator(HarmonicRotator&&) noexcept = delete;
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
}

#include "SpecialFunctionsImpl/Bessel.h"
#include "SpecialFunctionsImpl/Gamma.h"
#include "SpecialFunctionsImpl/Legendre.h"
#include "SpecialFunctionsImpl/Hermite.h"
#include "SpecialFunctionsImpl/IncompBeta.h"
