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
#pragma once

#include <QtWidgets/QWidget>
#include <QtWidgets/QVBoxLayout>
#include <QtDataVisualization/Q3DSurface>
#include <QtDataVisualization/QSurfaceDataProxy>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/LValueMatrix.h"

namespace Physica::Gui {
    class Plot3D : public QWidget {
        QVBoxLayout* vLayout;
    public:
        Plot3D(QWidget* parent = nullptr);
        /* Operations */
        template<class MatrixType>
        Q3DSurface& surf(const Core::LValueMatrix<MatrixType>& x,
                         const Core::LValueMatrix<MatrixType>& y,
                         const Core::LValueMatrix<MatrixType>& z);
    };

    template<class MatrixType>
    Q3DSurface& Plot3D::surf(const Core::LValueMatrix<MatrixType>& x,
                             const Core::LValueMatrix<MatrixType>& y,
                             const Core::LValueMatrix<MatrixType>& z) {
        assert(x.getColumn() == y.getColumn() && x.getRow() == y.getRow());
        assert(y.getColumn() == z.getColumn() && y.getRow() == z.getRow());
        auto* series = new QSurface3DSeries(new QSurfaceDataProxy());

        auto* dataArray = new QSurfaceDataArray;
        dataArray->reserve(z.getRow());
        for (size_t i = 0 ; i < z.getRow() ; ++i) {
            auto* dataRow = new QSurfaceDataRow(z.getColumn());
            for (size_t j = 0; j < z.getColumn(); j++)
                (*dataRow)[j].setPosition(QVector3D(float(x(i, j)),
                                                    float(z(i, j)),
                                                    float(y(i, j))));
            *dataArray << dataRow;
        }
        series->dataProxy()->resetArray(dataArray);

        auto* surface = new Q3DSurface();
        surface->addSeries(series);
        surface->axisX()->setLabelAutoRotation(30);
        surface->axisY()->setLabelAutoRotation(90);
        surface->axisZ()->setLabelAutoRotation(30);

        auto* widget = QWidget::createWindowContainer(surface, this);
        widget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        layout()->addWidget(widget);
        return *surface;
    }
}
