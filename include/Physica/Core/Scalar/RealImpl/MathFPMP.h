/*
 * Copyright 2024 Weibo He.
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

#include "FloatMP.h"

namespace Physica::Core {
    template<>
    __host__ __device__ Real<FloatMP> abs(const Real<FloatMP>& s) noexcept;

    template<>
    __host__ __device__ Real<FloatMP> sqrt(const Real<FloatMP>& s) noexcept;

    template<>
    __host__ __device__ Real<FloatMP> ln(const Real<FloatMP>& s) noexcept;

    template<>
    __host__ __device__ Real<FloatMP> exp(const Real<FloatMP>& s) noexcept;

    //!Compute a ^ unit.
    inline Real<FloatMP> powWord(const Real<FloatMP>& a, MPUnit unit) {
        Real<FloatMP> result(a);
        const auto lastUnitBits = std::countl_zero(unit);
        for(unsigned int j = 0; j < MPUnitWidth - lastUnitBits; ++j) {
            result = square(result);
            if((unit & 1U) != 0)
                result *= a;
            unit >>= 1U;
        }
        return result;
    }
    //!Compute a ^ unit, the highest bit of unit must be set.
    inline Real<FloatMP> powFullWord(const Real<FloatMP>& a, MPUnit unit) {
        Real<FloatMP> result(a);
        for(int j = 0; j < 64; ++j) {
            result = square(result);
            if((unit & 1U) != 0)
                result *= a;
            unit >>= 1U;
        }
        return result;
    }
    /*!
     * Calculate a^n.
     *
     * Reference: MaTHmu Project Group.计算机代数系统的数学原理[M].Beijing: TsingHua University Press, 2009:45
     */
    inline Real<FloatMP> powScalar(const Real<FloatMP>& a, const Real<FloatMP>& n) {
        const auto size = n.getSize();
        Real<FloatMP> result(a);
        if(n.getLength() < 0)
            result = reciprocal(a);

        for(int i = 0; i < size - 1; ++i)
            result = powFullWord(result, n[i]);
        result = powWord(result, n[size - 1]);
        return result;
    }

    template<>
    inline Real<FloatMP> pow(const Real<FloatMP>& s1, const Real<FloatMP>& s2) noexcept {
        return s1.isInteger() ? powScalar(s1, s2) : exp(ln(s1) * s2);
    }

    template<>
    Real<FloatMP> factorial(const Real<FloatMP>& s) noexcept;

    template<>
    __host__ __device__ Real<FloatMP> cos(const Real<FloatMP>& s) noexcept;

    template<>
    __host__ __device__ Real<FloatMP> sin(const Real<FloatMP>& s) noexcept;

    template<>
    Real<FloatMP> arccos(const Real<FloatMP>& s) noexcept;

    template<>
    Real<FloatMP> arcsin(const Real<FloatMP>& s) noexcept;

    template<>
    Real<FloatMP> arctan(const Real<FloatMP>& s) noexcept;

    template<>
    Real<FloatMP> cosh(const Real<FloatMP>& s) noexcept;

    template<>
    Real<FloatMP> sinh(const Real<FloatMP>& s) noexcept;

    template<>
    Real<FloatMP> tanh(const Real<FloatMP>& s) noexcept;

    template<>
    Real<FloatMP> sech(const Real<FloatMP>& s) noexcept;

    template<>
    Real<FloatMP> csch(const Real<FloatMP>& s) noexcept;

    template<>
    Real<FloatMP> coth(const Real<FloatMP>& s) noexcept;

    template<>
    Real<FloatMP> arccosh(const Real<FloatMP>& s) noexcept;

    template<>
    Real<FloatMP> arcsinh(const Real<FloatMP>& s) noexcept;

    template<>
    Real<FloatMP> arctanh(const Real<FloatMP>& s) noexcept;

    template<>
    Real<FloatMP> arccoth(const Real<FloatMP>& s) noexcept;

    template<>
    __host__ __device__ Real<FloatMP> floor(const Real<FloatMP>& s) noexcept;
}
