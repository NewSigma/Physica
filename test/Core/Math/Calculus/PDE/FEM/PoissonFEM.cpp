/*
 * Copyright 2022-2025 Weibo He.
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
#include "Physica/Core/Math/Calculus/PDE/FEM/Element/Rectangle1.h"
#include "Physica/Core/Math/Calculus/PDE/FEM/Element/Triangle1.h"
#include "Physica/Core/Math/Calculus/PDE/FEM/PoissonModel.h"
#include "Test.h"

using namespace Physica;
using T = float64;

constexpr double width = 2;
constexpr double height = 1;
constexpr double error = 1E-4;

namespace {
    T theory_solution(Vector2D<T> p) {
        constexpr double factor = -8 * width * width / (M_PI * M_PI * M_PI);
        const T b_2 = height * 0.5;
        const T rep_a = 1 / width;
        T result = 0;
        unsigned int k = 1;
        T temp = std::numeric_limits<T>::max();
        while (abs(temp) > error) {
            const T phase = k * M_PI;
            temp = cosh(phase * (p[1] - b_2) * rep_a) / (cosh(phase * b_2 * rep_a) * (k * k * k)) * sin(phase * p[0] * rep_a);
            k += 2;
            result += temp;
        }
        return result * factor - p[0] * (p[0] - width);
    }

    template<class ElementType>
    struct ElementIntegratorPacker {
        static T run(std::invocable<Vector2D<T>> auto fn) {
            return ElementType::gauss_integral(std::move(fn));
        }
    };
}

int main() {
    {
        using ElementType = Rectangle1<T>;
        auto mesh = ElementType::rectangle({0, 0}, {width, height}, 21, 21);
        mesh.addDirichletBoundary([](Vector2D<T> p) { return scalarNear(p[0], T(0), 1E-5)
                                                          || scalarNear(p[0], T(width), 1E-5)
                                                          || scalarNear(p[1], T(0), 1E-5)
                                                          || scalarNear(p[1], T(height), 1E-5); },
                                  [](Vector2D<T>) { return T(0); });

        auto func = [](Vector2D<T>) { return T(-2); };
        PoissonModel model(std::move(mesh), func);
        model.solve<ElementIntegratorPacker<ElementType>>();

        const VectorND<T> xs = VectorND<T>::linspace(0, width * 0.9, 6);
        const VectorND<T> ys = VectorND<T>::linspace(0, height * 0.9, 4);
        for (auto x : xs) {
            for (auto y : ys) {
                const T theory = theory_solution({x, y});
                const T simulation = model({x, y});
                expect(scalarNear(theory, simulation, 1E-2));
            }
        }
    }
    {
        using ElementType = Triangle1<T>;
        auto mesh = ElementType::rectangle({0, 0}, {width, height}, 20, 20);
        mesh.addDirichletBoundary([](Vector2D<T> p) { return scalarNear(p[0], T(0), 1E-5)
                                                          || scalarNear(p[0], T(width), 1E-5)
                                                          || scalarNear(p[1], T(0), 1E-5)
                                                          || scalarNear(p[1], T(height), 1E-5); },
                                  [](Vector2D<T>) { return T(0); });

        auto func = [](Vector2D<T>) { return T(-2); };
        PoissonModel model(std::move(mesh), func);
        model.solve<ElementIntegratorPacker<ElementType>>();

        const VectorND<T> xs = VectorND<T>::linspace(0, width * 0.9, 6);
        const VectorND<T> ys = VectorND<T>::linspace(0, height * 0.9, 4);
        for (auto x : xs) {
            for (auto y : ys) {
                const T theory = theory_solution({x, y});
                const T simulation = model({x, y});
                expect(scalarNear(theory, simulation, 1E-2));
            }
        }
    }
    return 0;
}
