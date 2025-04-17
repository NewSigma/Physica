/*
 * Copyright 2024-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/Core/Exception/BadConvergenceException.h"
#include "../RValueMatrix.h"

namespace Physica {
    template<Scalar T, size_t RowAtCompile, size_t ColAtCompile>
    class SVD;

    template<class Derived>
    auto RValueMatrix<Derived>::norm1() const -> Tr {
        Tr maxNorm1 = 0;
        for (size_t c = 0; c < getCol(); ++c) {
            const auto v = col(c);
            Tr temp = abs(v).sum();
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
    auto RValueMatrix<Derived>::norm1_power(unsigned int maxIteration) const -> Tr {
        assert(getRow() == getCol() && "[Error]: norm1_power only applies to square matrix");
        assert(maxIteration > 0 && "[Error]: Invalid max iteration");
        using VectorType = VectorND<ScalarType>;
        const Derived& m = Base::getDerived();
        const size_t length = getRow();
        const Trv factor = reciprocal(Trv(length));

        VectorND<Tr> x(length, factor);
        VectorType y(length);
        VectorType z(length);
        unsigned int iteration = 0;
        size_t lastIndex = 0;
        size_t index = 0;
        while (iteration < maxIteration) {
            y.template operator=<decltype(m * x), Executor>(m * x);

            using MatVecDot = decltype(m.hermite() * unit(y));
            z.template operator=<MatVecDot, Executor>(m.hermite() * unit(y));
            const Trv criteria = iteration == 0 ? (z.reals().sum().value() * factor) : z[index].real().value();
            const bool isConverged = z.values().normInf() <= criteria * Trv(std::numeric_limits<T>::epsilon() + 1);
            const bool isCycling = iteration > 0 && (lastIndex == index); // Avoid cycling because unit(0) is implemented as 1
            if (isConverged || isCycling) {
                if constexpr (isComplex)
                    return y.norm1();
                else {
                    VectorType v = VectorType::linspace(1, 2, length) + ScalarType(1);
                    for (size_t i = 0; i < length; i += 2)
                        v[i] = -v[i];
                    x = m * v;
                    const Tr pick = x.norm1() / v.norm1();
                    return std::max(y.norm1(), pick);
                }
            }

            if (iteration == 0)
                x = Tr(0);
            else
                x[index] = Tr(0);

            size_t nextIndex = 0;
            Tr maxAbs = 0;
            for (size_t i = 0; i < length; ++i) {
                const Tr temp = abs(z.calc(i).real());
                if (temp > maxAbs) {
                    nextIndex = i;
                    maxAbs = temp;
                }
            }
            if (lastIndex != nextIndex)
                lastIndex = index;
            index = nextIndex;
            x[index] = Trv(1);
            iteration += 1;
        }
        throw BadConvergenceException("[Error]: norm1_power failed to converge");
    }

    template<class Derived>
    auto RValueMatrix<Derived>::normF() const -> Tr {
        return sqrt(squaredNorms().sum());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::normInf() const -> Tr {
        Trv result = std::numeric_limits<T>::lowest();
        for (size_t i = 0; i < getRow(); ++i)
            result = std::max(abs(row(i)).sum(), result);
        return result;
    }

    template<class Derived>
    auto RValueMatrix<Derived>::cond2() const -> T {
        SVD<T, RowAtCompile, ColAtCompile> svd(getRow(), getCol());
        svd.compute(Base::getDerived());
        const auto& s = svd.getSingulars();
        return s.max() / s.min();
    }
}
