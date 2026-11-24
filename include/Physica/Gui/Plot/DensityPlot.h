/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2024 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <QtWidgets/QWidget>
#include "Physica/Core/Math/Geometry/Point.h"

namespace Physica {
    /**
     * DensityPlot provide a approach to show 2D data.
     */
    class PHYSICA_API DensityPlot : public QWidget {
    private:
        static constexpr int GlobalAlpha = 200;

        QImage dataImage;
        QImage colorImage;
        QRgb colorList[256]{};
    public:
        explicit DensityPlot(int width, int height, QWidget* parent = nullptr);
        ~DensityPlot() noexcept override = default;
        /* Operations */
        void appendPoint(int x, int y, int radius, unsigned char alpha);
    protected:
        void paintEvent(QPaintEvent* event) override;
    private:
        void initColorList();
        void drawColorImage();
    };
}
