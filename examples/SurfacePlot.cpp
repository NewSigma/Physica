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
#include <QApplication>
#include "Physica/Gui/Plot/Plot3D.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

using namespace Physica::Core;
using namespace Physica::Gui;
using ScalarType = Scalar<Double, false>;

int main(int argc, char **argv) {
    QApplication app(argc, argv);

    Vector<ScalarType> v = Vector<ScalarType>::linspace(-5, 5, 50);
    auto grid = DenseMatrix<ScalarType>::meshgrid(v, v);
    DenseMatrix<ScalarType> z = hadamard(square(grid.first) + square(grid.second), exp(-square(grid.first) - square(grid.second) + ScalarType::One()));

    Plot3D* plot3d = new Plot3D();

    auto& surf = plot3d->surf(grid.first, grid.second, z);
    surf.activeTheme()->setType(Q3DTheme::ThemePrimaryColors);

    {
        QLinearGradient gr;
        gr.setColorAt(0.0, Qt::blue);
        gr.setColorAt(0.35, Qt::cyan);
        gr.setColorAt(0.5, Qt::green);
        gr.setColorAt(0.65, Qt::yellow);
        gr.setColorAt(1.0, Qt::red);

        surf.seriesList().at(0)->setBaseGradient(gr);
        surf.seriesList().at(0)->setColorStyle(Q3DTheme::ColorStyleRangeGradient);
    }
    plot3d->show();

    return QApplication::exec();
}
