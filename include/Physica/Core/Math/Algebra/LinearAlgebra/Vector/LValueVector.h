/*
 * Copyright 2021 WeiBo He.
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

#include "RValueVector.h"

namespace Physica::Core {
    /**
     * \class LValueVector is base class of vectors that can be assigned to \class LValueVector
     * and other vectors can be assigned to this class.
     * In other words, you can take the address of elements in the vector.
     */
    template<class Derived>
    class LValueVector : public RValueVector<Derived> {
        using Base = RValueVector<Derived>;
    public:
        using typename Base::ScalarType;
    public:
        [[nodiscard]] ScalarType& operator[](size_t index) { return Base::getDerived()[index]; }
        [[nodiscard]] const ScalarType& operator[](size_t index) const { return Base::getDerived()[index]; }
    };
}
