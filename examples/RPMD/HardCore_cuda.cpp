/*
 * Copyright 2023 WeiBo He.
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
#include <fstream>
#include <gperftools/profiler.h>
#include <QApplication>
#include <QtCharts/QValueAxis>
#include "Physica/Gui/Plot/Plot.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"
#include "Physica/Core/Math/Statistics/ProbabilityDistributionFunction.h"
#include "Physica/Core/Math/Statistics/NumCharacter.h"
#include "Physica/Utils/Random.h"

using namespace Physica::Core;
using namespace Physica::Gui;
using ScalarType = Scalar<Double, false>;
using VectorType = Vector<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;

void run(unsigned int sys, std::mt19937& gen, const ProbabilityDistributionFunction<ScalarType>& originPdf, MatrixType& record);

int main(int argc, char** argv) {
    const ProbabilityDistributionFunction<ScalarType> pdf(-5, 5, 200);
    MatrixType record(200, 8);
    //ProfilerStart("profiler.dat");
    SequentialExecutor::parallel_for([&record, &pdf](unsigned int sys) {
        std::mt19937::result_type seed;
        Physica::Utils::Random::rdrand(seed);
        std::mt19937 gen(seed);

        run(sys, gen, pdf, record);
    }, 8, 1).wait();
    return 0;

    VectorType mean(record.getRow()), devia(record.getRow());
    for (size_t i = 0; i < mean.getLength(); ++i) {
        mean[i] = Physica::Core::mean(record.row(i));
        devia[i] = Physica::Core::deviation(record.row(i));
    }

    std::ofstream fout("data");
    fout << mean << devia;

    QApplication app(argc, argv);
    QFont font;
    Plot* plot = new Plot();
    auto& chart = *plot->chart();
    chart.legend()->setVisible(true);
    chart.legend()->setAlignment(Qt::AlignTop);
    chart.legend()->setMarkerShape(QLegend::MarkerShapeFromSeries);
    {
        double minX = double(pdf.getFromPoint());
        double maxX = double(pdf.getToPoint());
        constexpr double minY = 0;
        constexpr double maxY = 0.4;
        constexpr double deltaX = 2;
        constexpr double deltaY = 0.1;
        QValueAxis* axisX = new QValueAxis();
        font = axisX->labelsFont();
        font.setPointSize(15);
        axisX->setTickAnchor(0);
        axisX->setTickInterval(deltaX);
        axisX->setTickType(QValueAxis::TicksDynamic);
        axisX->setMinorGridLineVisible(false);
        axisX->setLinePenColor(Qt::black);
        axisX->setGridLineVisible(false);
        axisX->setLabelsFont(font);
        axisX->setRange(minX, maxX);
        axisX->setTitleText("v");
        axisX->setTitleFont(font);
        QValueAxis* axisY = new QValueAxis();
        axisY->setTickAnchor(0);
        axisY->setTickInterval(deltaY);
        axisY->setTickType(QValueAxis::TicksDynamic);
        axisY->setMinorGridLineVisible(false);
        axisY->setMinorTickCount(4);
        axisY->setLinePenColor(Qt::black);
        axisY->setGridLineVisible(false);
        axisY->setMinorGridLineVisible(false);
        axisY->setLabelsFont(font);
        axisY->setRange(minY, maxY);
        axisY->setTitleText("P(v)");
        axisY->setTitleFont(font);
        QValueAxis* axisTop = new QValueAxis();
        axisTop->setTickAnchor(0);
        axisTop->setTickInterval(deltaX);
        axisTop->setTickType(QValueAxis::TicksDynamic);
        axisTop->setLabelsVisible(false);
        axisTop->setGridLineVisible(false);
        axisTop->setRange(minX, maxX);
        axisTop->setLinePenColor(Qt::black);
        QValueAxis* axisRight = new QValueAxis();
        axisRight->setTickAnchor(0);
        axisRight->setTickInterval(deltaY);
        axisRight->setTickType(QValueAxis::TicksDynamic);
        axisRight->setLabelsVisible(false);
        axisRight->setGridLineVisible(false);
        axisRight->setMinorGridLineVisible(false);
        axisRight->setMinorTickCount(4);
        axisRight->setRange(minY, maxY);
        axisRight->setLinePenColor(Qt::black);

        chart.addAxis(axisX, Qt::AlignBottom);
        chart.addAxis(axisY, Qt::AlignLeft);
        chart.addAxis(axisTop, Qt::AlignTop);
        chart.addAxis(axisRight, Qt::AlignRight);

        VectorType x = pdf.makePotision();
        {
            auto& area = plot->area_center(x, mean, devia);
            area.attachAxis(axisX);
            area.attachAxis(axisY);

            auto& spline = plot->line(x, mean);
            spline.setColor(area.color());
            spline.attachAxis(axisX);
            spline.attachAxis(axisY);

            auto color = area.color();
            color.setAlpha(75);
            area.setColor(color);
            spline.setName("Simulation");
        }
        {
            ScalarType mass = M_SQRT2;
            ScalarType temperatureT = 2;
            VectorType y = sqrt(mass / (temperatureT * (2 * M_PI))) * exp(-square(x) * (mass / (temperatureT * 2.0)));
            auto& spline = plot->line(x, y);
            auto pen = spline.pen();
            pen.setColor(Qt::blue);
            pen.setStyle(Qt::DashLine);
            spline.setPen(pen);
            spline.attachAxis(axisX);
            spline.attachAxis(axisY);
            spline.setName("Boltzmann");
        }
    }
    plot->show();
    return QApplication::exec();
}
