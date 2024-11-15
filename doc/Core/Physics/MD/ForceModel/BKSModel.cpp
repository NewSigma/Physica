/*
 * Copyright 2024 Weibo He.
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
#include <QApplication>
#include <Physica/Core/Physics/MD/ForceModel/Ewald/Ewald.h>
#include <Physica/Core/Physics/MD/ForceModel/BKSModel.h>
#include <Physica/Gui/Plot/Plot.h>

using namespace Physica::Core;
using namespace Physica::Gui;
using ScalarType = Real<Double>;
using VectorType = VectorND<ScalarType>;
using ForceModel = BKSModel<ScalarType, Ewald<ScalarType, RSpaceEwald<ScalarType>>, false>;
constexpr double A = ForceModel::A_SiO;
constexpr double b = ForceModel::b_SiO;
constexpr double c = ForceModel::c_SiO;

ScalarType pot_functor(ScalarType r) {
    const ScalarType r2 = square(r);
    const ScalarType r4 = square(r2);
    return ScalarType(A) * exp(ScalarType(-b) * r) - ScalarType(c) / (r2 * r4);
}

ScalarType force_functor(ScalarType r) {
    const ScalarType r2 = square(r);
    const ScalarType r4 = square(r2);
    return ScalarType(A * b) * exp(ScalarType(-b) * r) - ScalarType(6 * c) / (r * r2 * r4);
}

int main(int argc, char** argv) {
    const VectorType x = VectorType::linspace(0.01, 5, 500);
    VectorType pot(x.getLength());
    VectorType force(x.getLength());
    for (size_t i = 0; i < x.getLength(); ++i) {
        pot[i] = pot_functor(x[i]);
        force[i] = force_functor(x[i]);
    }

    QApplication app(argc, argv);
    {
        Plot* plot = new Plot(1, 5, -2, 1, 1, 1);
        plot->getChart()->legend()->setVisible(false);
        auto* axisX = plot->getAxisX();
        auto* axisLeft = plot->getAxisY();
        auto* axisRight = plot->getAxisRight();
        axisX->setTitleText("r/Bohr");
        axisX->setLabelFormat("%d");
        axisLeft->setTitleText("Energy/Hartree");
        axisLeft->setLabelFormat("%d");
        axisRight->setTitleText("Force/a.u.");
        axisRight->setLabelFormat("%d");
        axisRight->setLabelsVisible(true);
        {
            auto& line = plot->line(VectorType{0, 5}, VectorType{0, 0});
            auto pen = line.pen();
            pen.setColor(Qt::black);
            pen.setStyle(Qt::DashLine);
            line.setPen(pen);
        }
        {
            const auto color = Qt::blue;
            auto pen = axisLeft->linePen();
            pen.setColor(color);
            axisLeft->setLinePen(pen);
            axisLeft->setLabelsColor(color);

            auto& line = plot->line(x, pot);
            line.setColor(color);
        }
        {
            const auto color = Qt::red;
            auto pen = axisRight->linePen();
            pen.setColor(color);
            axisRight->setLinePen(pen);
            axisRight->setLabelsColor(color);

            auto& line = plot->line(x, force);
            line.setColor(color);
            line.attachAxis(axisRight);
        }
        plot->show();
    }
    return QApplication::exec();
}
