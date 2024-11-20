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
#include <QApplication>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Core;
using namespace Physica::Gui;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;

int main(int argc, char** argv) {
    const auto betas = VectorType::linspace(0, 4, 41);
    VectorType tpq_hphi, fulldiag_hphi, tpq, fulldiag;
    H5File h5f("TPQ_FD.h5", H5File::ReadOnly);
    tpq_hphi.read(h5f, "TPQ_HPhi");
    fulldiag_hphi.read(h5f, "FullDiag_HPhi");
    tpq.read(h5f, "TPQ");
    fulldiag.read(h5f, "FullDiag");

    QApplication app(argc, argv);
    Plot* plot = new Plot(0, 4, 0.5, 1.05, 1, 0.2);
    plot->getLegend().setAlignment(Qt::AlignRight);
    auto* axisX = plot->getAxisX();
    auto* axisY = plot->getAxisY();
    axisX->setTitleText("&beta;");
    axisX->setLabelFormat("%d");
    axisY->setTitleText("&rho;");
    axisY->setLabelFormat("%.1f");

    auto& scatter_tpq = plot->scatter(betas, tpq_hphi);
    scatter_tpq.setName("TPQ(HPhi)");
    scatter_tpq.setMarkerSize(10);
    auto pen1 = scatter_tpq.pen();
    pen1.setColor(scatter_tpq.color());
    scatter_tpq.setPen(pen1);
    scatter_tpq.setColor(Qt::transparent);

    auto& scatter_fd = plot->scatter(betas, fulldiag_hphi);
    scatter_fd.setName("FullDiag(HPhi)");
    scatter_fd.setMarkerSize(10);
    auto pen2 = scatter_fd.pen();
    pen2.setColor(scatter_fd.color());
    scatter_fd.setPen(pen2);
    scatter_fd.setColor(Qt::transparent);

    auto& line_tpq = plot->line(betas, tpq);
    line_tpq.setName("TPQ");
    line_tpq.setColor(pen1.color());

    auto& line_fd = plot->line(betas, fulldiag);
    line_fd.setName("FullDiag");
    auto pen = line_fd.pen();
    pen.setColor(pen2.color());
    pen.setStyle(Qt::DashLine);
    QList<qreal> dashes;
    dashes << 20 << 20;
    pen.setDashPattern(dashes);
    line_fd.setPen(pen);

    plot->show();
    return QApplication::exec();
}
