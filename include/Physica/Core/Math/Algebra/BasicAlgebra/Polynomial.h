/*
 * Copyright 2019-2021 WeiBo He.
 *
 * This file is part of Physica.

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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/RValueVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/EigenSolver.h"

namespace Physica::Core {
    /**
     * Class Polynomial provides operations around polynomials.
     *
     * A polynomial is:
     * y(x) = a0 + a1 x + a2 x ^ 2 + ... + (an - 1) x ^ (n - 1) + 1 x^n
     */
    template<class ScalarType, size_t Power>
    class Polynomial {
        using VectorType = Vector<ScalarType, Power>;
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
        template<class AnyScalar>
        typename Internal::BinaryScalarOpReturnType<ScalarType, AnyScalar>::Type operator()(const ScalarBase<AnyScalar>& x) const;
        /* Getters */
        [[nodiscard]] size_t getPower() const noexcept { return coeffs.getLength(); }
        [[nodiscard]] const VectorType& getCoeffVector() const noexcept { return coeffs; }
    };

    template<class ScalarType, size_t Power>
    template<class AnyScalar>
    typename Internal::BinaryScalarOpReturnType<ScalarType, AnyScalar>::Type
    Polynomial<ScalarType, Power>::operator()(const ScalarBase<AnyScalar>& x) const {
        using ResultType = typename Internal::BinaryScalarOpReturnType<ScalarType, AnyScalar>::Type;
        if (coeffs.empty())
            return ResultType::Zero();
        ResultType result = ResultType(coeffs[0]);
        AnyScalar temp = x.getDerived();
        const auto length = coeffs.getLength();
        for (size_t i = 1; i < length; ++i) {
            result += temp * coeffs[i];
            temp *= x.getDerived();
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
    template<class ScalarType, size_t Power>
    Vector<ComplexScalar<ScalarType>, Power>
    polyRoot(const Polynomial<ScalarType, Power>& poly) {
        using MatrixType = DenseMatrix<ScalarType, DenseMatrixOption::Column | DenseMatrixOption::Vector, Power, Power>;

        const size_t power = poly.getPower();
        MatrixType companion = MatrixType::Zeros(power);
        for (size_t i = 0; i < power - 1; ++i)
            companion(i + 1, i) = ScalarType::One();
        auto col = companion.col(power - 1);
        col = -poly.getCoeffVector();

        EigenSolver<MatrixType> solver(companion, false);
        return solver.getEigenvalues();
    }
}
