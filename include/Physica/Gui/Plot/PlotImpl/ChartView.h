/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
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
