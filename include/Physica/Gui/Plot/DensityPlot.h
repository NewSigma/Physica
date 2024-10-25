/*
 * Copyright 2020-2024 Weibo He.
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
#include <Physica/Core/Math/Geometry/Point.h>

namespace Physica::Gui {
    /*!
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
