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

#include "Element.h"

namespace Physica::Core {
    template<class Scalar>
    class Square4P : public Element<2, Scalar> {
    public:
        ~Square4P() = default;
        /* Operations */
        Scalar shapePartialS1(ShapeIndex shapeIndex, const Poing<2>& p) override final;
        Scalar shapePartialS2(ShapeIndex shapeIndex, const Poing<2>& p) override final;
        Matrix jacobi(const Poing<2>& p) override final;
    protected:
        Square4P() : Element(4) {}
    };
}

#include "Square4PImpl.h"
