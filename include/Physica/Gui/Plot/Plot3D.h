/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
        QSurface3DSeries& surf(const Matrix auto& x, const Matrix auto& y, const Matrix auto& z);
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

    QSurface3DSeries& Plot3D::surf(const Matrix auto& x, const Matrix auto& y, const Matrix auto& z) {
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
