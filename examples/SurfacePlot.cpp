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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Transform/FFT.h"
#include "Physica/Gui/Plot/Plot3D.h"

using namespace Physica::Core;
using namespace Physica::Gui;
using namespace Physica::Utils;
using RealType = Scalar<Double, false>;
using ComplexType = ComplexScalar<RealType>;

int main(int argc, char** argv) {
    const size_t N1 = 100;
    const size_t N2 = 100;
    const double deltaX = 0.01;
    const double deltaY = 0.01;

    Vector<ComplexType> data(N1 * N2);
    {
        size_t index = 0;
        for (size_t i = 0; i < N1; ++i) {
            for (size_t j = 0; j < N2; ++j) {
                data[index++] = ComplexType(RealType(std::sin(2 * M_PI * 10 * i * deltaX) + 2 * std::cos(2 * M_PI * 5 * j * deltaY)), RealType::Zero());
            }
        }
    }
    FFT<ComplexType, 2> fft(data, {N1, N2}, {deltaX, deltaY});
    fft.transform();
    auto components = fft.getComponents();

    DenseMatrix<RealType> f(N1 / 2 * 2 + 1, N2 / 2 * 2 + 1);
    {
        size_t index = 0;
        for (size_t i = 0; i < N1 / 2 * 2 + 1; ++i) {
            for (size_t j = 0; j < N2 / 2 * 2 + 1; ++j) {
                f(i, j) = components[index].norm();
                index++;
            }
        }
    }

    QApplication app(argc, argv);

    Vector<RealType> v = Vector<RealType>::linspace(-1 / (2 * deltaX), 1 / (2 * deltaX), N1 / 2 * 2 + 1);
    auto grid = DenseMatrix<RealType>::meshgrid(v, v);

    Plot3D* plot3d = new Plot3D();

    auto& surf = plot3d->surf(grid.first, grid.second, f);
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

    return app.exec();
}
