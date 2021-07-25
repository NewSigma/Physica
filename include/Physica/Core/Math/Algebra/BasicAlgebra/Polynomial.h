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

#include "Physica/Core/MultiPrecition/Scalar.h"
#include "Physica/Core/Utils/Container/Array/Array.h"

namespace Physica::Core {
    /*!
     * Class Polynomial provides operations around polynomials.
     *
     * A polynomial is:
     * y(x) = a0 + a1 x + a2 x ^ 2 + ... + an x ^ n
     */
    template<class ScalarType>
    class Polynomial {
    public:
        //data stores the coeficients of the polynomial.
        Array<ScalarType> data;
    public:
        Polynomial() = default;
        Polynomial(const Polynomial& p) = default;
        Polynomial(Polynomial&& p) noexcept : data(std::move(p.data)) {}
        ~Polynomial() = default;
        /* Operators */
        Polynomial& operator=(const Polynomial& p) = default;
        Polynomial& operator=(Polynomial&& p) noexcept { data = std::move(p.data); return *this; }
        ScalarType operator()(const ScalarType& s) const;
    };

    /*!
     * Returns the value of this polynomial when x = s.
     */
    template<class ScalarType>
    ScalarType Polynomial<ScalarType>::operator()(const ScalarType& s) const {
        if(data.empty())
            return ScalarType::Zero();
        ScalarType result(data[0]);
        ScalarType temp(s);
        const auto length = data.getLength();
        for(size_t i = 1; i < length; ++i) {
            result += temp * data[i];
            temp *= s;
        }
        return result;
    }
}
