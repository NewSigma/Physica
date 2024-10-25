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

namespace Physica::Core {
    template<>
    Scalar<FloatMP> abs(const Scalar<FloatMP>& s) noexcept;

    template<>
    Scalar<FloatMP> sqrt(const Scalar<FloatMP>& s) noexcept;

    template<>
    Scalar<FloatMP> ln(const Scalar<FloatMP>& s) noexcept;

    template<>
    Scalar<FloatMP> exp(const Scalar<FloatMP>& s) noexcept;

    //!Compute a ^ unit.
    inline Scalar<FloatMP> powWord(const Scalar<FloatMP>& a, MPUnit unit) {
        Scalar<FloatMP> result(a);
        const auto lastUnitBits = countLeadingZeros(unit);
        for(unsigned int j = 0; j < MPUnitWidth - lastUnitBits; ++j) {
            result = square(result);
            if((unit & 1U) != 0)
                result *= a;
            unit >>= 1U;
        }
        return result;
    }
    //!Compute a ^ unit, the highest bit of unit must be set.
    inline Scalar<FloatMP> powFullWord(const Scalar<FloatMP>& a, MPUnit unit) {
        Scalar<FloatMP> result(a);
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
    inline Scalar<FloatMP> powScalar(const Scalar<FloatMP>& a, const Scalar<FloatMP>& n) {
        const auto size = n.getSize();
        Scalar<FloatMP> result(a);
        if(n.getLength() < 0)
            result = reciprocal(a);

        for(int i = 0; i < size - 1; ++i)
            result = powFullWord(result, n[i]);
        result = powWord(result, n[size - 1]);
        return result;
    }

    template<>
    inline Scalar<FloatMP> pow(const Scalar<FloatMP>& s1, const Scalar<FloatMP>& s2) noexcept {
        return s1.isInteger() ? powScalar(s1, s2) : exp(ln(s1) * s2);
    }

    template<>
    Scalar<FloatMP> factorial(const Scalar<FloatMP>& s) noexcept;

    template<>
    Scalar<FloatMP> cos(const Scalar<FloatMP>& s) noexcept;

    template<>
    Scalar<FloatMP> sin(const Scalar<FloatMP>& s) noexcept;

    template<>
    Scalar<FloatMP> arccos(const Scalar<FloatMP>& s) noexcept;

    template<>
    Scalar<FloatMP> arcsin(const Scalar<FloatMP>& s) noexcept;

    template<>
    Scalar<FloatMP> arctan(const Scalar<FloatMP>& s) noexcept;

    template<>
    Scalar<FloatMP> cosh(const Scalar<FloatMP>& s) noexcept;

    template<>
    Scalar<FloatMP> sinh(const Scalar<FloatMP>& s) noexcept;

    template<>
    Scalar<FloatMP> tanh(const Scalar<FloatMP>& s) noexcept;

    template<>
    Scalar<FloatMP> sech(const Scalar<FloatMP>& s) noexcept;

    template<>
    Scalar<FloatMP> csch(const Scalar<FloatMP>& s) noexcept;

    template<>
    Scalar<FloatMP> coth(const Scalar<FloatMP>& s) noexcept;

    template<>
    Scalar<FloatMP> arccosh(const Scalar<FloatMP>& s) noexcept;

    template<>
    Scalar<FloatMP> arcsinh(const Scalar<FloatMP>& s) noexcept;

    template<>
    Scalar<FloatMP> arctanh(const Scalar<FloatMP>& s) noexcept;

    template<>
    Scalar<FloatMP> arccoth(const Scalar<FloatMP>& s) noexcept;

    template<>
    Scalar<FloatMP> floor(const Scalar<FloatMP>& s) noexcept;
}
