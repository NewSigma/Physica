/*
 * Copyright 2021-2023 Weibo He.
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

namespace Physica {
    class PHYSICA_API Plot3D : public QWidget {
        QVBoxLayout* vLayout;
    protected:
        Q3DSurface* surface;
    public:
        Plot3D(QWidget* parent = nullptr);
        /* Operations */
        template<Matrix T>
        QSurface3DSeries& surf(const T& x,
                               const T& y,
                               const T& z);
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

    template<Matrix T>
    QSurface3DSeries& Plot3D::surf(const T& x,
                                   const T& y,
                                   const T& z) {
        assert(x.getCol() == y.getCol() && x.getRow() == y.getRow());
        assert(y.getCol() == z.getCol() && y.getRow() == z.getRow());
        auto* series = new QSurface3DSeries(new QSurfaceDataProxy());

        auto* dataArray = new QSurfaceDataArray;
        dataArray->reserve(z.getRow());
        for (size_t i = 0 ; i < z.getRow() ; ++i) {
            auto* dataRow = new QSurfaceDataRow(z.getCol());
            for (size_t j = 0; j < z.getCol(); j++)
                (*dataRow)[j].setPosition(QVector3D(float(x.calc(i, j)),
                                                    float(z.calc(i, j)),
                                                    float(y.calc(i, j))));
            *dataArray << dataRow;
        }
        series->dataProxy()->resetArray(dataArray);
        surface->addSeries(series);
        return *series;
    }
}
