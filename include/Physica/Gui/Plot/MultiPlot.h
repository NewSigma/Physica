/*
 * Copyright 2024-2025 Weibo He.
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

#include <QGridLayout>
#include "Plot.h"

namespace Physica {
    class PHYSICA_API MultiPlot : public QWidget {
        using This = MultiPlot;
        using Base = QWidget;

        QGridLayout* layout;
    public:
        MultiPlot(int row, int col, QWidget* parent = nullptr);
        MultiPlot(int row, int col, double minX, double maxX, double minY, double maxY, double deltaX, double deltaY, QWidget* parent = nullptr);
        MultiPlot(const This&) = delete;
        MultiPlot(This&&) noexcept = delete;
        ~MultiPlot() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        Plot& operator[](int row, int col);
        /* Operations */
        void forPlot(std::invocable<Plot&> auto fn) noexcept;

        void toSvg(const char* path, int resolution = 0);
        /* Getters */
        [[nodiscard]] int getRow() const noexcept { return layout->rowCount(); }
        [[nodiscard]] int getCol() const noexcept { return layout->columnCount(); }
        /* Setters */
        void setMinX(double value) noexcept { forPlot([=](Plot& p) noexcept { p.setMinX(value); }); }
        void setMaxX(double value) noexcept { forPlot([=](Plot& p) noexcept { p.setMaxX(value); }); }
        void setRangeX(double minX, double maxX) { forPlot([=](Plot& p) noexcept { p.setRangeX(minX, maxX); }); }
        void setMinY(double value) noexcept { forPlot([=](Plot& p) noexcept { p.setMinY(value); }); }
        void setMaxY(double value) noexcept { forPlot([=](Plot& p) noexcept { p.setMaxY(value); }); }
        void setRangeY(double minY, double maxY) { forPlot([=](Plot& p) noexcept { p.setRangeY(minY, maxY); }); }
        void setDeltaX(double value) noexcept { forPlot([=](Plot& p) noexcept { p.setDeltaX(value); }); }
        void setDeltaY(double value) noexcept { forPlot([=](Plot& p) noexcept { p.setDeltaY(value); }); }
        inline void setBox(double minX, double maxX, double minY, double maxY, double deltaX, double deltaY);

        void setFont(QFont f) { forPlot([&](Plot& p) { p.setFont(f); }); }
        void setFontSize(int size) { forPlot([&](Plot& p) { p.setFontSize(size); }); }
    };

    void MultiPlot::forPlot(std::invocable<Plot&> auto fn) noexcept {
        for (int r = 0; r < getRow(); ++r)
            for (int c = 0; c < getCol(); ++c)
                fn(operator[](r, c));
    }

    inline void MultiPlot::setBox(double minX, double maxX, double minY, double maxY, double deltaX, double deltaY) {
        forPlot([=](Plot& p) noexcept { p.setBox(minX, maxX, minY, maxY, deltaX, deltaY); });
    }
}
