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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/UnitVector.h"
#include "Physica/Core/Exception/BadConvergenceException.h"
#include "../RValueMatrix.h"

namespace Physica {
    template<Scalar T, size_t RowAtCompile, size_t ColAtCompile>
    class SVD;

    template<class Derived>
    auto RValueMatrix<Derived>::norm1() const -> Tr {
        Tr maxNorm1 = 0;
        for (size_t c = 0; c < getCol(); ++c) {
            const auto v = Base::getDerived().col(c);
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
    template<ExecutePolicy P>
    auto RValueMatrix<Derived>::norm1_power(unsigned int maxIteration) const -> Tr {
        assert(isSquare() && "[Error]: norm1_power only applies to square matrix");
        assert(maxIteration > 0 && "[Error]: Invalid max iteration");
        const size_t length = getRow();
        if (length == 1) [[unlikely]]
            return abs(calc(0, 0));

        using Buffer = DenseVector<T, RowAtCompile>;
        const Trv factor = reciprocal(Trv(length));
        Buffer y(length);
        Buffer z(length);
        size_t lastIndex = 0;
        size_t index = 0;
        for (unsigned int iteration = 0; iteration < maxIteration; ++iteration) {
            const Derived& m = Base::getDerived();
            if (iteration == 0) {
                z = factor;
                (m * z).template assign<P>(y);
            }
            else
                (m * UnitVector<Trv>(index, length)).template assign<P>(y);
            (m.hermite() * unit(y)).template assign<P>(z);

            const Trv criteria = iteration == 0 ? (z.reals().sum().value() * factor) : z[index].real().value();
            const bool isConverged = z.values().normInf() <= criteria * Trv(std::numeric_limits<T>::epsilon() + 1);
            const bool isCycling = iteration > 0 && (lastIndex == index); // Avoid cycling because unit(0) is implemented as 1
            if (isConverged || isCycling) {
                const Tr normY = y.norm1();
                if constexpr (isComplex())
                    return normY;
                else {
                    y.linspace(1, 2);
                    for (size_t i = 1; i < length; i += 2)
                        y[i] = -y[i];
                    z = m * y;
                    return std::max(normY, z.norm1() / y.norm1());
                }
            }

            const size_t nextIndex = [&z, length]() noexcept {
                size_t result{};
                Tr zmax = 0;
                for (size_t i = 0; i < length; ++i) {
                    const Tr temp = isComplex() ? z[i].squaredNorm() : abs(z[i]);
                    if (temp > zmax) {
                        result = i;
                        zmax = temp;
                    }
                }
                assert(zmax.isPositive() && "[Error]: Unexpected empty vector");
                return result;
            }();

            if (lastIndex != nextIndex)
                lastIndex = index;
            index = nextIndex;
        }
        throw BadConvergenceException("[Error]: norm1_power failed to converge");
    }

    template<class Derived>
    auto RValueMatrix<Derived>::normF() const -> CoDiff<Tr> {
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
