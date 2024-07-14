/*
 * Copyright 2024 WeiBo He.
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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h>
#include <Physica/Core/Exception/BadConvergenceException.h>

namespace Physica::Core {
    template<class Derived>
    typename RValueMatrix<Derived>::RealType RValueMatrix<Derived>::norm1() const {
        RealType maxNorm1 = 0;
        for (size_t c = 0; c < getColumn(); ++c) {
            const auto v = col(c);
            RealType temp = abs(v).sum();
            if (temp > maxNorm1)
                maxNorm1 = temp;
        }
        return maxNorm1;
    }
    /**
     * Estimate the norm-1 using power method
     * 
     * Reference:
     * [1] SIAM J. Sci. Stat. Comput. 11, 804–809 (1990); https://doi.org/10.1137/0911047
     */
    template<class Derived>
    template<class Executor>
    typename RValueMatrix<Derived>::RealType RValueMatrix<Derived>::norm1_power(unsigned int maxIteration) const {
        assert(getRow() == getColumn() && "[Error]: norm1_power only applies to square matrix");
        assert(maxIteration > 0 && "[Error]: Invalid max iteration");
        using VectorType = Vector<ScalarType>;
        const Derived& m = Base::getDerived();
        const size_t length = getRow();
        const RealType factor = reciprocal(RealType(length));

        Vector<RealType> x(length, factor);
        VectorType y(length);
        VectorType z(length);
        unsigned int iteration = 0;
        size_t index = 0;
        while (iteration < maxIteration) {
            y.template operator=<decltype(m * x), Executor>(m * x);

            using ExprType = decltype(m.hermite() * unit(y));
            z.template operator=<ExprType, Executor>(m.hermite() * unit(y));
            const RealType criteria = iteration == 0 ? (toRealVector(z).sum() * factor) : z[index].getReal();
            if (z.normInf() <= criteria) {
                if constexpr (isComplex)
                    return y.norm1();
                else {
                    VectorType v = VectorType::linspace(0, 1, length) + ScalarType(1);
                    for (size_t i = 0; i < length; i += 2)
                        v[i] = -v[i];
                    x = m * v;
                    const RealType pick = x.norm1() / v.norm1();
                    return std::max(y.norm1(), pick);
                }
            }

            if (iteration == 0)
                x = RealType(0);
            else
                x[index] = RealType(0);

            RealType maxAbs = 0;
            for (size_t i = 0; i < length; ++i) {
                const RealType temp = abs(z.calc(i).getReal());
                if (temp > maxAbs) {
                    index = i;
                    maxAbs = temp;
                }
            }
            x[index] = RealType(1);
            iteration += 1;
        }
        throw BadConvergenceException("[Error]: norm1_power failed to converge");
    }
}
