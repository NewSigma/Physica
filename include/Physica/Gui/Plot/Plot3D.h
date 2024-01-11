/*
 * Copyright 2021-2023 WeiBo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/LValueMatrix.h"

namespace Physica::Gui {
    class PHYSICA_API Plot3D : public QWidget {
        QVBoxLayout* vLayout;
    protected:
        Q3DSurface* surface;
    public:
        Plot3D(QWidget* parent = nullptr);
        /* Operations */
        template<class MatrixType>
        QSurface3DSeries& surf(const Core::LValueMatrix<MatrixType>& x,
                               const Core::LValueMatrix<MatrixType>& y,
                               const Core::LValueMatrix<MatrixType>& z);
        /* Getters */
        [[nodiscard]] QValue3DAxis* getAxisX() const noexcept { return surface->axisX(); }
        [[nodiscard]] QValue3DAxis* getAxisY() const noexcept { return surface->axisZ(); }
        [[nodiscard]] QValue3DAxis* getAxisZ() const noexcept { return surface->axisY(); }
        [[nodiscard]] float getMinZ() const noexcept { return getAxisZ()->min(); }
        [[nodiscard]] float getMaxZ() const noexcept { return getAxisZ()->max(); }
        /* Setters */
        void setMinZ(float value) { getAxisZ()->setMin(value); }
        void setMaxZ(float value) { getAxisZ()->setMax(value); }
        /* Static member */
        [[nodiscard]] static QLinearGradient makeDefaultGrad();
    };

    template<class MatrixType>
    QSurface3DSeries& Plot3D::surf(const Core::LValueMatrix<MatrixType>& x,
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
        surface->addSeries(series);
        return *series;
    }
}
