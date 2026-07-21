/*
 * Copyright 2019-2026 Weibo He.
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
module;

#include <QtCharts/QValueAxis>
#include <QtCharts/QLineSeries>
#include <QtCharts/QSplineSeries>
#include <QtCharts/QScatterSeries>
#include <QtCharts/QAreaSeries>
#include <QtCharts/QBoxPlotSeries>
#include <QtCharts/QLegend>
#include "Physica/Core/Scalar/Real.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

export module Physica.Gui.Plot;

export import Physica.Gui.ChartView;
export import Physica.Gui.ContourSeries;

export namespace Physica {
    class PHYSICA_API Plot : public ChartView {
        using Base = ChartView;

        QValueAxis* axisX;
        QValueAxis* axisY;
        QValueAxis* axisTop;
        QValueAxis* axisRight;
    public:
        Plot(QWidget* parent = nullptr);
        Plot(float64 minX, float64 maxX, float64 minY, float64 maxY, float64 deltaX, float64 deltaY, QWidget* parent = nullptr);
        Plot(const Plot&) = delete;
        Plot(Plot&&) noexcept = delete;
        ~Plot() = default;
        /* Operators */
        Plot& operator=(const Plot&) = delete;
        Plot& operator=(Plot&&) noexcept = delete;
        /* Operations */
        QLineSeries& line(const Vector auto& y);
        QLineSeries& line(const Vector auto& x, const Vector auto& y);
        QSplineSeries& spline(const Vector auto& y);
        QSplineSeries& spline(const Vector auto& x, const Vector auto& y);
        QScatterSeries& scatter(const Vector auto& y);
        QScatterSeries& scatter(const Vector auto& x, const Vector auto& y);
        QAreaSeries& hist(const Vector auto& data, size_t binCount, bool density = false);
        QAreaSeries& area_boundary(const Vector auto& x,
                                   const Vector auto& lower,
                                   const Vector auto& upper);
        QAreaSeries& area_center(const Vector auto& x,
                                 const Vector auto& center,
                                 const Vector auto& deviation);
        template<Matrix M>
        ContourSeries<M>& contour(const M& x, const M& y, const M& z, Array<double> levels);
        template<Vector V>
        QBoxPlotSeries& boxWhisker(const V& x, const Array<V>& data);
        QBoxPlotSeries& errorBar(const Vector auto& mean, const Vector auto& deviation);
        QBoxPlotSeries& errorBar(const Vector auto& x, const Vector auto& mean, const Vector auto& deviation);
        QScatterSeries& label(double x, double y, QString text);

        [[nodiscard]] QImage toImage(int width = 640, int height = 480);
        void toSvg(const char* path, int width = 640, int height = 480, int resolution = 0);
        /* Getters */
        [[nodiscard]] QLegend& getLegend() const noexcept { return *getChart()->legend(); }
        [[nodiscard]] QValueAxis* getAxisX() const noexcept { return axisX; }
        [[nodiscard]] QValueAxis* getAxisY() const noexcept { return axisY; }
        [[nodiscard]] QValueAxis* getAxisTop() const noexcept { return axisTop; }
        [[nodiscard]] QValueAxis* getAxisRight() const noexcept { return axisRight; }
        [[nodiscard]] double getMinX() const noexcept { return axisX->min(); }
        [[nodiscard]] double getMaxX() const noexcept { return axisX->max(); }
        [[nodiscard]] double getMinY() const noexcept { return axisY->min(); }
        [[nodiscard]] double getMaxY() const noexcept { return axisY->max(); }
        /* Setters */
        void setAxisX(QValueAxis* axis);
        void setAxisY(QValueAxis* axis);
        void setAxisTop(QValueAxis* axis);
        void setAxisRight(QValueAxis* axis);
        void setMinX(double value) noexcept;
        void setMaxX(double value) noexcept;
        void setRangeX(double minX, double maxX);
        void setMinY(double value) noexcept;
        void setMaxY(double value) noexcept;
        void setRangeY(double minY, double maxY);
        void setDeltaX(double value) noexcept;
        void setDeltaY(double value) noexcept;
        void setBox(float64 minX, float64 maxX, float64 minY, float64 maxY, float64 deltaX, float64 deltaY);

        void setTickDirection(QAbstractAxis::TickDirection d);
        void setFont(const QFont& font);
        void setFontSize(int size);
    private:
        QBoxSet* setFromVector(const Vector auto& v);
        double findMedian(const Vector auto& sorted_v, size_t from, size_t to);
    };

    QLineSeries& Plot::line(const Vector auto& y) {
        return line(VectorND<float64>::linspace(0, y.getLength() - 1, y.getLength()), y);
    }

    QLineSeries& Plot::line(const Vector auto& x, const Vector auto& y) {
        assert(x.getLength() == y.getLength());
        auto* series = new QLineSeries();
        for (size_t i = 0; i < x.getLength(); ++i)
            *series << QPointF(double(x.calc(i)), double(y.calc(i)));
        getChart()->addSeries(series);
        series->attachAxis(axisX);
        series->attachAxis(axisY);

        update();
        return *series;
    }

    QSplineSeries& Plot::spline(const Vector auto& y) {
        return spline(VectorND<float64>::linspace(0, y.getLength() - 1, y.getLength()), y);
    }

    QSplineSeries& Plot::spline(const Vector auto& x, const Vector auto& y) {
        assert(x.getLength() == y.getLength());
        auto* series = new QSplineSeries();
        for (size_t i = 0; i < x.getLength(); ++i)
            *series << QPointF(double(x.calc(i)), double(y.calc(i)));
        getChart()->addSeries(series);
        series->attachAxis(axisX);
        series->attachAxis(axisY);

        update();
        return *series;
    }

    QScatterSeries& Plot::scatter(const Vector auto& y) {
        return scatter(VectorND<float64>::linspace(0, y.getLength() - 1, y.getLength()), y);
    }

    QScatterSeries& Plot::scatter(const Vector auto& x, const Vector auto& y) {
        assert(x.getLength() == y.getLength());
        auto* series = new QScatterSeries();
        for (size_t i = 0; i < x.getLength(); ++i)
            *series << QPointF(double(x.calc(i)), double(y.calc(i)));
        getChart()->addSeries(series);
        series->attachAxis(axisX);
        series->attachAxis(axisY);

        update();
        return *series;
    }

    QAreaSeries& Plot::hist(const Vector auto& data, size_t binCount, bool density) {
        double binWidth{};
        double min{};
        const size_t length = data.getLength();
        /* Get binWidth and min */ {
            auto minimum = data.calc(0);
            auto maximum = data.calc(0);
            for (size_t i = 1; i < length; ++i) {
                auto temp = data.calc(i);
                if (temp < minimum)
                    minimum = std::move(temp);
                else if (temp > maximum)
                    maximum = std::move(temp);
            }
            assert(maximum >= minimum);
            min = double(minimum);
            binWidth = double(maximum - minimum + (binCount - 1)) / double(binCount);
            if (binWidth == 0)
                binWidth = 1;
        }

        Array<unsigned int> arr(binCount + 1, 0);
        const double binCountPerUnit = 1 / binWidth;
        for (size_t i = 0; i < length; ++i) {
            const auto binIndex = size_t((double(data.calc(i)) - min) * binCountPerUnit);
            arr[binIndex]++;
        }

        auto* upper_series = new QLineSeries();
        const double initial_x = min;
        double current_x = initial_x;
        if (density) {
            const double density_factor = 1 / (binWidth * length);
            for (size_t i = 0; i < binCount; ++i) {
                const double y = arr[i] * density_factor;
                *upper_series << QPointF(current_x, y);
                current_x += binWidth;
                *upper_series << QPointF(current_x, y);
            }
        }
        else {
            for (size_t i = 0; i < binCount; ++i) {
                const double y = double(arr[i]);
                *upper_series << QPointF(current_x, y);
                current_x += binWidth;
                *upper_series << QPointF(current_x, y);
            }
        }
        auto* lower_series = new QLineSeries();
        *lower_series << QPointF(initial_x, 0) << QPointF(current_x, 0);

        auto* series = new QAreaSeries(upper_series, lower_series);

        getChart()->addSeries(series);
        series->attachAxis(axisX);
        series->attachAxis(axisY);

        update();
        return *series;
    }

    QAreaSeries& Plot::area_boundary(const Vector auto& x, const Vector auto& lower, const Vector auto& upper) {
        assert(x.getLength() == lower.getLength() && x.getLength() == upper.getLength());

        auto* upper_series = new QLineSeries();
        auto* lower_series = new QLineSeries();
        for (size_t i = 0; i < x.getLength(); ++i) {
            const double x_i = double(x.calc(i));
            *lower_series << QPointF(x_i, double(lower.calc(i)));
            *upper_series << QPointF(x_i, double(upper.calc(i)));
        }
        auto* series = new QAreaSeries(upper_series, lower_series);
        getChart()->addSeries(series);
        series->attachAxis(axisX);
        series->attachAxis(axisY);
        update();
        return *series;
    }

    QAreaSeries& Plot::area_center(const Vector auto& x, const Vector auto& center, const Vector auto& deviation) {
        return area_boundary(x, center - deviation, center + deviation);
    }

    template<Matrix M>
    ContourSeries<M>& Plot::contour(const M& x, const M& y, const M& z, Array<double> levels) {
        auto* series = new ContourSeries<M>(x, y, z, std::move(levels));
        series->attachTo(*getChart());
        series->attachAxis(axisX);
        series->attachAxis(axisY);
        update();
        return *series;
    }

    template<Vector V>
    QBoxPlotSeries& Plot::boxWhisker(const V& x, const Array<V>& data) {
        assert(x.getLength() == data.getLength());
        auto* series = new QBoxPlotSeries(QBoxPlotSeries::Numeric);
        for (size_t i = 0; i < x.getLength(); ++i) {
            auto* set = setFromVector(data[i]);
            set->setX(double(std::move(x.calc(i))));
            series->append(set);
        }
        getChart()->addSeries(series);
        series->attachAxis(axisX);
        series->attachAxis(axisY);
        update();
        return *series;
    }

    QBoxPlotSeries& Plot::errorBar(const Vector auto& mean, const Vector auto& deviation) {
        using namespace Physica;
        using ScalarType = float64;
        const auto x = VectorND<ScalarType>::linspace(0, mean.getLength() - 1, mean.getLength());
        return errorBar(x, mean, deviation);
    }

    QBoxPlotSeries& Plot::errorBar(const Vector auto& x, const Vector auto& mean, const Vector auto& deviation) {
        assert(x.getLength() == mean.getLength() && x.getLength() == deviation.getLength());
        auto* series = new QBoxPlotSeries(QBoxPlotSeries::Numeric);
        for (size_t i = 0; i < x.getLength(); ++i) {
            if (deviation.calc(i).isNegative() || !deviation.calc(i).isFinite())
                continue;
            auto* set = new QBoxSet();
            const auto mean_i = double(mean.calc(i));
            const auto devia_i = double(deviation.calc(i));
            set->setValue(QBoxSet::LowerExtreme, mean_i - devia_i);
            set->setValue(QBoxSet::UpperExtreme, mean_i + devia_i);
            set->setValue(QBoxSet::Median, mean_i);
            set->setValue(QBoxSet::LowerQuartile, mean_i);
            set->setValue(QBoxSet::UpperQuartile, mean_i);
            set->setX(double(x.calc(i)));
            series->append(set);
        }
        getChart()->addSeries(series);
        series->attachAxis(axisX);
        series->attachAxis(axisY);
        update();
        return *series;
    }

    QBoxSet* Plot::setFromVector(const Vector auto& v) {
        VectorND<float64> buffer = v;
        std::sort(buffer.begin(), buffer.end());
        auto* result = new QBoxSet();
        const size_t count = v.getLength();
        result->setValue(QBoxSet::LowerExtreme, double(buffer.front()));
        result->setValue(QBoxSet::UpperExtreme, double(buffer.back()));
        result->setValue(QBoxSet::Median, findMedian(buffer, 0, count));
        result->setValue(QBoxSet::LowerQuartile, findMedian(buffer, 0, count / 2));
        result->setValue(QBoxSet::UpperQuartile, findMedian(buffer, (count / 2) + (count % 2), count));
        return result;
    }

    double Plot::findMedian(const Vector auto& sorted_v, size_t from, size_t to) {
        size_t count = to - from;
        size_t count2 = count / 2;
        if (count % 2)
            return double(sorted_v[count2 + from]);

        auto right = sorted_v[count2 + from];
        auto left = sorted_v[count2 - 1 + from];
        return double((right + left) * 0.5);
    }
}
