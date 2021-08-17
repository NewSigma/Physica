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

namespace Physica::Core {
    template<class ScalarType>
    ComplexScalar<ScalarType> square(const ComplexScalar<ScalarType>& c);

    template<class ScalarType>
    inline ComplexScalar<ScalarType> reciprocal(const ComplexScalar<ScalarType>& c);

    template<class ScalarType>
    ComplexScalar<ScalarType> sqrt(const ComplexScalar<ScalarType>& c);

    template<class ScalarType>
    ComplexScalar<ScalarType> ln(const ComplexScalar<ScalarType>& c);

    template<class ScalarType>
    ComplexScalar<ScalarType> exp(const ComplexScalar<ScalarType>& c);

    template<class ScalarType>
    ComplexScalar<ScalarType> cos(const ComplexScalar<ScalarType>& c);

    template<class ScalarType>
    ComplexScalar<ScalarType> sin(const ComplexScalar<ScalarType>& c);

    template<class ScalarType>
    ComplexScalar<ScalarType> tan(const ComplexScalar<ScalarType>& c);

    template<class ScalarType>
    ComplexScalar<ScalarType> sec(const ComplexScalar<ScalarType>& c);

    template<class ScalarType>
    ComplexScalar<ScalarType> csc(const ComplexScalar<ScalarType>& c);

    template<class ScalarType>
    ComplexScalar<ScalarType> cot(const ComplexScalar<ScalarType>& c);

    template<class ScalarType>
    ComplexScalar<ScalarType> cosh(const ComplexScalar<ScalarType>& c);

    template<class ScalarType>
    ComplexScalar<ScalarType> sinh(const ComplexScalar<ScalarType>& c);

    template<class ScalarType>
    ComplexScalar<ScalarType> tanh(const ComplexScalar<ScalarType>& c);

    template<class ScalarType>
    ComplexScalar<ScalarType> sech(const ComplexScalar<ScalarType>& c);

    template<class ScalarType>
    ComplexScalar<ScalarType> csch(const ComplexScalar<ScalarType>& c);

    template<class ScalarType>
    ComplexScalar<ScalarType> coth(const ComplexScalar<ScalarType>& c);
}

#include "FunctionImpl/CElementaryImpl.h"
