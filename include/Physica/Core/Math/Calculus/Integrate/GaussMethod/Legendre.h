/*
 * Copyright 2022 WeiBo He.
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
    template<class ScalarType, unsigned int Dim, unsigned int Order>
    class GaussLegendre;

    template<class ScalarType, unsigned int Dim>
    class GaussLegendre<ScalarType, Dim, 1> {
        using VectorType = Vector<ScalarType, Dim>;
        constexpr static unsigned int Factor = 1U << Dim;
    public:
        template<class Functor>
        static ScalarType run(Functor func);
    };
    /**
     * \param func
     * ScalarType Functor(VectorType)
     */
    template<class ScalarType, unsigned int Dim>
    template<class Functor>
    ScalarType GaussLegendre<ScalarType, Dim, 1>::run(Functor func) {
        return func(VectorType(Dim, 0)) * Factor;
    }
}
