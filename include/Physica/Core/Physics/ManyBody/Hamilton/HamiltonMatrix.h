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

#include <iostream>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/RValueMatrix.h"

namespace Physica::Core {
    template<class Derived>
    class HamiltonMatrix : public RValueMatrix<Derived> {
        using This = HamiltonMatrix<Derived>;
        using Base = RValueMatrix<Derived>;
        using ReprType = typename Traits<Derived>::ReprType;
        constexpr static unsigned int NumSite = ReprType::StateType::NumSite;
    public:
        ~HamiltonMatrix() = default;
        /* Operations */
        [[nodiscard]] const This& hermite() const noexcept { return *this; }
        /* Getters */
        [[nodiscard]] size_t getNumState() const noexcept { return Base::getDerived().getNumState(); }
        [[nodiscard]] size_t getRow() const noexcept { return getNumState(); }
        [[nodiscard]] size_t getCol() const noexcept { return getNumState(); }
    protected:
        HamiltonMatrix() = default;
        HamiltonMatrix(const This&) = default;
        HamiltonMatrix(This&&) noexcept = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
    };
}

namespace Physica {
    template<class Derived>
    class Traits<Core::HamiltonMatrix<Derived>> {
    public:
        constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::Element;
        constexpr static size_t RowAtCompile = Dynamic;
        constexpr static size_t ColAtCompile = Dynamic;
        constexpr static size_t SizeAtCompile = Dynamic;
    };
}
