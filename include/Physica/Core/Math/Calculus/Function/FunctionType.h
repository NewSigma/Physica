/*
 * Copyright 2020 WeiBo He.
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
#ifndef PHYSICA_FUNCTIONTYPE_H
#define PHYSICA_FUNCTIONTYPE_H

namespace Physica::Core {
    enum FunctionType {
        Value,
        //functions
        Add,
        Sub,
        Mul,
        Div,
        Square,
        Reciprocal,
        Sqrt,
        Factorial,
        Ln,
        Log,
        Exp,
        Cos,
        Sin,
        Tan,
        Sec,
        Csc,
        Cot,
        ArcCos,
        ArcSin,
        ArcTan,
        ArcSec,
        ArcCsc,
        ArcCot,
        Cosh,
        Sinh,
        Tanh,
        Sech,
        Csch,
        Coth,
        ArcCosh,
        ArcSinh,
        ArcTanh,
        ArcSech,
        ArcCsch,
        ArcCoth,
    };
}

#endif
