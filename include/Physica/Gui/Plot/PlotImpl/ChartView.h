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
#pragma once

#include <QGraphicsView>
#include <QChart>
#include <QRubberBand>
#include "Physica/Macro.h"

namespace Physica {
    class PHYSICA_API ChartView : public QGraphicsView {
        using Base = QGraphicsView;

        QChart* chart;
        QRubberBand& rubberBand;
        QPoint rubberBandOrigin;
    public:
        explicit ChartView(QWidget* parent = nullptr) : ChartView(nullptr, parent) {}
        explicit ChartView(QChart* chart_, QWidget* parent = nullptr);
        ChartView(const ChartView&) = delete;
        ChartView(ChartView&&) noexcept = delete;
        ~ChartView() = default;
        /* Operators */
        ChartView& operator=(const ChartView&) = delete;
        ChartView& operator=(ChartView&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] QChart* getChart() const { return chart; }
        /* Setters */
        void setChart(QChart* chart_);
        void setRubberBandEnabled(bool isEnabled) { rubberBand.setEnabled(isEnabled); }
    protected:
        void resizeEvent(QResizeEvent* event) override;
        void mousePressEvent(QMouseEvent* event) override;
        void mouseMoveEvent(QMouseEvent* event) override;
        void mouseReleaseEvent(QMouseEvent* event) override;
    private:
        /* Operations */
        void resize();
    };
}
