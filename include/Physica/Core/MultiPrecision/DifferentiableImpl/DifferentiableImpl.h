/*
 * Copyright 2023 WeiBo He.
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
    inline Differentiable<ScalarType> Differentiable<ScalarType>::operator+(const Differentiable<ScalarType>& other) const {
        return {value + other.value, tangent + other.tangent};
    }

    template<class ScalarType>
    inline Differentiable<ScalarType> Differentiable<ScalarType>::operator-(const Differentiable<ScalarType>& other) const {
        return {value - other.value, tangent - other.tangent};
    }

    template<class ScalarType>
    inline Differentiable<ScalarType> Differentiable<ScalarType>::operator*(const Differentiable<ScalarType>& other) const {
        return {value * other.value, tangent * other.value + value * other.tangent};
    }

    template<class ScalarType>
    inline Differentiable<ScalarType> Differentiable<ScalarType>::operator/(const Differentiable<ScalarType>& other) const {
        const ScalarType rep = reciprocal(other.value);
        return {value * rep, (tangent * other.value - value * other.tangent) * square(rep)}
    }

    template<class ScalarType>
    inline Differentiable<ScalarType> Differentiable<ScalarType>::operator-() const {
        return {-value, -tangent};
    }
}
