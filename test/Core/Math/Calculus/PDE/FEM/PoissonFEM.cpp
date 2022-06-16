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
#include "Physica/Utils/TestHelper.h"
#include "Physica/Core/Math/Calculus/PDE/FEM/PoissonModel.h"
#include "Physica/Core/Math/Calculus/PDE/FEM/Element/Rectangle1.h"

using namespace Physica::Core;
using ScalarType = Scalar<Double, false>;
using VectorType = Vector<ScalarType, 2>;
using ElementType = Rectangle1<ScalarType>;

constexpr double width = 2;
constexpr double height = 1;
constexpr double error = 1E-4;

ScalarType theory_solution(VectorType p) {
    constexpr double factor = -8 * width * width / (M_PI * M_PI * M_PI);
    const ScalarType b_2 = height * 0.5;
    const ScalarType rep_a = 1 / width;
    ScalarType result = 0;
    unsigned int k = 1;
    ScalarType temp = std::numeric_limits<ScalarType>::max();
    while (abs(temp) > error) {
        const ScalarType phase = k * M_PI;
        temp = cosh(phase * (p[1] -b_2) * rep_a) / (cosh(phase * b_2 * rep_a) * (k * k * k)) * sin(phase * p[0] * rep_a);
        k += 2;
        result += temp;
    }
    return result * factor - p[0] * (p[0] - width);
}

struct ElementIntegratorPacker {
    template<class Functor>
    static ScalarType run(Functor func) {
        return ElementType::integral(std::move(func));
    }
};

int main() {
    auto mesh = ElementType::rectangle({0, 0}, {width, height}, 21, 21);
    mesh.addDirichletBoundary([](VectorType p) { return scalarNear(p[0], ScalarType::Zero(), 1E-5)
                                                      || scalarNear(p[0], ScalarType(width), 1E-5)
                                                      || scalarNear(p[1], ScalarType::Zero(), 1E-5)
                                                      || scalarNear(p[1], ScalarType(height), 1E-5); },
                              []([[maybe_unused]] VectorType p) { return ScalarType(0); });

    PoissonModel model(std::move(mesh));
    auto func = []([[maybe_unused]] VectorType p) { return ScalarType(-2); };
    model.solve<decltype(func), ElementIntegratorPacker>(func);

    const Vector<ScalarType> xs = Vector<ScalarType>::linspace(0, width * 0.9, 6);
    const Vector<ScalarType> ys = Vector<ScalarType>::linspace(0, height * 0.9, 4);
    ScalarType max_relative_error = 0;
    for (auto x : xs) {
        for (auto y : ys) {
            const ScalarType theory = theory_solution({x, y});
            const ScalarType simulation = model({x, y});
            if (!scalarNear(theory, simulation, 1E-2))
                return 1;
        }
    }
    return 0;
}
