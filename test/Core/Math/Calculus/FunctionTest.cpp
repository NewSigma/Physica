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
#include "Physica/Core/Math/Calculus/Function/VectorFunction/VectorFunction.h"
#include <functional>

using namespace Physica::Core;

int main() {
    TreeFunction<Double, false> func(1, 2);
    /* Initialize func */ {
        func.setConstant(2, 0);
        func.setConstant(3.14159, 1);
        auto tree = func.getConstantNode(1) * func.getVariableNode(0) / func.getConstantNode(0) * sin(func.getConstantNode(1) * func.getVariableNode(0));
        func.setTree(std::make_shared<decltype(tree)>(std::move(tree)));
    }
    std::function<Scalar<Double, false>(Scalar<Double, false>)> a(func);
    double d1 = a(0.5).getTrivial();
    double d2 = func(0.5).getTrivial();
    return (d1 != d2);
}