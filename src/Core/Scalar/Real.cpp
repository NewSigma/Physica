/*
 * Copyright 2026 Weibo He.
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
#include <print>
#include "Physica/Core/Scalar/Real.h"

using namespace Physica;

namespace {
    template<FloatPrec Prec>
    Real<Prec> stripSignificandImpl(Real<Prec> x) noexcept {
        assert(!x.isSubNormal() && "[Error]: Input must be a normal positive value");
        int exp = 0;
        std::frexp(x.toMachine(), &exp);
        return Real<Prec>(std::ldexp(Real<Prec>(1).toMachine(), exp - 1));
    }
}

auto Real<Float32>::stripSignificand() const noexcept -> Real {
    return stripSignificandImpl(*this);
}

void Real<Float32>::dump() const noexcept {
    std::println("{}", *this);
}

auto Real<Float64>::stripSignificand() const noexcept -> Real {
    return stripSignificandImpl(*this);
}

void Real<Float64>::dump() const noexcept {
    std::println("{}", *this);
}
