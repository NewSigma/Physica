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
#include <iostream>
#include "Physica/Core/Math/Transform/DFT.h"
#include <QtWidgets/QApplication>
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Core;
using namespace Physica::Gui;
using ScalarType = Scalar<Double, false>;

ScalarType func(ScalarType x) {
    return sin(ScalarType(2 * M_PI * 3) * x) + sin(ScalarType(6 * M_PI) * x) * 5;
}

int main(int argc, char** argv) {
    const size_t N = 100;
    const double t_max = 1;
    const Vector<ScalarType> x = Vector<ScalarType>::linspace(ScalarType::Zero(), ScalarType(t_max), N);
    Vector<ScalarType> data(N);
    for (size_t i = 0; i < N; ++i)
        data[i] = func(x[i]);

    DFT<ScalarType> dft(data, ScalarType(t_max / N));
    dft.transform();
    Vector<ScalarType> f = Vector<ScalarType>::linspace(ScalarType(N / t_max / -2), ScalarType(N / t_max / 2), N);
    Vector<ScalarType> data1 = toImagVector(dft.getComponents());

    QApplication app(argc, argv);
    Plot* plot = new Plot();
    auto& scatter = plot->scatter(f, data1);
    scatter.setMarkerSize(8);
    plot->chart()->createDefaultAxes();
    plot->show();
    return QApplication::exec();
}
