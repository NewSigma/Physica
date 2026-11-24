/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2019-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/EigenSolver.h"

namespace Physica {
    /**
     * Class Polynomial provides operations around polynomials.
     *
     * A polynomial is:
     * y(x) = a0 + a1 x + a2 x ^ 2 + ... + (an - 1) x ^ (n - 1) + 1 x^n
     */
    template<Scalar T, size_t Power>
    class Polynomial {
        using VectorType = DenseVector<T, Power>;
    private:
        VectorType coeffs;
    public:
        Polynomial() = default;
        Polynomial(const VectorType& coeffs_) : coeffs(coeffs_) {}
        Polynomial(const Polynomial& p) = default;
        Polynomial(Polynomial&& p) noexcept : coeffs(std::move(p.coeffs)) {}
        ~Polynomial() = default;
        /* Operators */
        Polynomial& operator=(const Polynomial& p) = default;
        Polynomial& operator=(Polynomial&& p) noexcept { coeffs = std::move(p.coeffs); return *this; }
        auto operator()(const auto& x) const;
        /* Getters */
        [[nodiscard]] size_t getPower() const noexcept { return coeffs.getLength(); }
        [[nodiscard]] const VectorType& getCoeffVector() const noexcept { return coeffs; }
    };

    template<Scalar T, size_t Power>
    auto Polynomial<T, Power>::operator()(const auto& x) const {
        using U = std::remove_cvref<decltype(x)>::type;
        using ResultType = Internal::BinaryScalarOpRtnTy<T, U>::Type;
        if (coeffs.empty())
            return ResultType(0);
        ResultType result = ResultType(coeffs[0]);
        U temp = x;
        const auto length = coeffs.getLength();
        for (size_t i = 1; i < length; ++i) {
            result += temp * coeffs[i];
            temp *= x;
        }
        result += temp;
        return result;
    }
    /**
     * Solve polynomial equation
     *
     * Reference:
     * [1] https://www.mathworks.com/help/matlab/ref/roots.html
     */
    template<Scalar T, size_t Power>
    DenseVector<Complex<T>, Power>
    polyRoot(const Polynomial<T, Power>& poly) {
        using MatrixType = DenseMatrix<T, MatrixMajor::Col, Power, Power>;

        const size_t power = poly.getPower();
        auto companion = MatrixType::zeros(power);
        for (size_t i = 0; i < power - 1; ++i)
            companion[i + 1, i] = T(1);
        auto col = companion.col(power - 1);
        col = -poly.getCoeffVector();

        EigenSolver<T, Power> solver(companion, false);
        return solver.getEigenvalues();
    }
}
