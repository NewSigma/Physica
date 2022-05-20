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

using namespace Physica::Core;
using ScalarType = Scalar<Double, false>;
using VectorType = Vector<ScalarType, 2>;

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

int main() {
    auto mesh = rectangle<ScalarType>({0, 0}, {width, height}, 20, 10);
    mesh.addDirichletBoundary([](VectorType p) { return scalarNear(p[0], ScalarType::Zero(), 1E-5)
                                                      || scalarNear(p[0], ScalarType(width), 1E-5)
                                                      || scalarNear(p[1], ScalarType::Zero(), 1E-5)
                                                      || scalarNear(p[1], ScalarType(height), 1E-5); },
                              []([[maybe_unused]] VectorType p) { return ScalarType(0); });

    PoissonModel model(std::move(mesh));
    auto func = []([[maybe_unused]] VectorType p) { return ScalarType(-2); };
    model.solve<decltype(func), typename decltype(model)::Integrator>(func);

    const Vector<ScalarType> xs = Vector<ScalarType>::linspace(0, width * 0.9, 6);
    const Vector<ScalarType> ys = Vector<ScalarType>::linspace(0, height * 0.9, 4);
    ScalarType rmes = 0;
    for (auto x : xs) {
        for (auto y : ys)
            rmes += square(theory_solution({x, y}) - model({x, y}));
    }
    rmes = sqrt(rmes / ScalarType(xs.getLength() * ys.getLength()));
    std::cout << rmes << std::endl;
    const bool isTrueAnswer = scalarNear(rmes, ScalarType::Zero(), error * 10);
    return !isTrueAnswer;
}
