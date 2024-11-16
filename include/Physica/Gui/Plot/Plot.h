/*
 * Copyright 2019-2024 Weibo He.
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

#include <QtCharts/QValueAxis>
#include <QtCharts/QLineSeries>
#include <QtCharts/QSplineSeries>
#include <QtCharts/QScatterSeries>
#include <QtCharts/QAreaSeries>
#include <QtCharts/QBoxPlotSeries>
#include <Physica/Core/MultiPrecision/Real.h>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h>
#include "PlotImpl/ChartView.h"
#include "ContourSeries.h"

namespace Physica::Gui {
    class PHYSICA_API Plot : public ChartView {
        using Base = ChartView;

        QValueAxis* axisX;
        QValueAxis* axisY;
        QValueAxis* axisTop;
        QValueAxis* axisRight;
    public:
        Plot(QWidget* parent = nullptr);
        Plot(double minX, double maxX, double minY, double maxY, double deltaX, double deltaY, QWidget* parent = nullptr);
        Plot(const Plot&) = delete;
        Plot(Plot&&) noexcept = delete;
        ~Plot() = default;
        /* Operators */
        Plot& operator=(const Plot&) = delete;
        Plot& operator=(Plot&&) noexcept = delete;
        /* Operations */
        template<class VectorType>
        QLineSeries& line(const Core::RValueVector<VectorType>& y);
        template<class VectorType1, class VectorType2>
        QLineSeries& line(const Core::RValueVector<VectorType1>& x, const Core::RValueVector<VectorType2>& y);
        template<class VectorType>
        QSplineSeries& spline(const Core::RValueVector<VectorType>& y);
        template<class VectorType1, class VectorType2>
        QSplineSeries& spline(const Core::RValueVector<VectorType1>& x, const Core::RValueVector<VectorType2>& y);
        template<class VectorType>
        QScatterSeries& scatter(const Core::RValueVector<VectorType>& y);
        template<class VectorType1, class VectorType2>
        QScatterSeries& scatter(const Core::RValueVector<VectorType1>& x, const Core::RValueVector<VectorType2>& y);
        template<class VectorType>
        QAreaSeries& hist(const Core::RValueVector<VectorType>& data, size_t binCount, bool density = false);
        template<class VectorType>
        QAreaSeries& area_boundary(const Core::RValueVector<VectorType>& x,
                                   const Core::RValueVector<VectorType>& lower,
                                   const Core::RValueVector<VectorType>& upper);
        template<class VectorType>
        QAreaSeries& area_center(const Core::RValueVector<VectorType>& x,
                                 const Core::RValueVector<VectorType>& center,
                                 const Core::RValueVector<VectorType>& deviation);
        template<class MatrixType>
        ContourSeries<MatrixType>& contour(const Core::LValueMatrix<MatrixType>& x,
                                           const Core::LValueMatrix<MatrixType>& y,
                                           const Core::LValueMatrix<MatrixType>& z,
                                           Core::Array<double> levels);
        template<class VectorType>
        QBoxPlotSeries& boxWhisker(const Core::RValueVector<VectorType>& x, const Core::Array<VectorType>& data);
        template<class VectorType>
        QBoxPlotSeries& errorBar(const Core::RValueVector<VectorType>& mean, const Core::RValueVector<VectorType>& deviation);
        template<class VectorType1, class VectorType2>
        QBoxPlotSeries& errorBar(const Core::RValueVector<VectorType1>& x, const Core::RValueVector<VectorType2>& mean, const Core::RValueVector<VectorType2>& deviation);
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
        inline void setAxisX(QValueAxis* axis);
        inline void setAxisY(QValueAxis* axis);
        inline void setAxisTop(QValueAxis* axis);
        inline void setAxisRight(QValueAxis* axis);
        inline void setMinX(double value) noexcept;
        inline void setMaxX(double value) noexcept;
        inline void setRangeX(double minX, double maxX);
        inline void setMinY(double value) noexcept;
        inline void setMaxY(double value) noexcept;
        inline void setRangeY(double minY, double maxY);
        inline void setDeltaX(double value) noexcept;
        inline void setDeltaY(double value) noexcept;
        void setBox(double minX, double maxX, double minY, double maxY, double deltaX, double deltaY);

        void setTickDirection(QAbstractAxis::TickDirection d);
        void setFont(QFont font);
        void setFontSize(int size);
    private:
        template<class VectorType>
        QBoxSet* setFromVector(const Core::LValueVector<VectorType>& v);
        template<class VectorType>
        double findMedian(const Core::LValueVector<VectorType>& sorted_v, size_t from, size_t to);
    };

    template<class VectorType>
    QLineSeries& Plot::line(const Core::RValueVector<VectorType>& y) {
        using Vector = Core::DenseVector<typename VectorType::ScalarType, VectorType::SizeAtCompile>;
        return line(Vector::linspace(0, y.getLength() - 1, y.getLength()), y);
    }

    template<class VectorType1, class VectorType2>
    QLineSeries& Plot::line(const Core::RValueVector<VectorType1>& x, const Core::RValueVector<VectorType2>& y) {
        assert(x.getLength() == y.getLength());
        QLineSeries* series = new QLineSeries();
        for (size_t i = 0; i < x.getLength(); ++i)
            *series << QPointF(double(x.calc(i)), double(y.calc(i)));
        getChart()->addSeries(series);
        series->attachAxis(axisX);
        series->attachAxis(axisY);

        update();
        return *series;
    }

    template<class VectorType>
    QSplineSeries& Plot::spline(const Core::RValueVector<VectorType>& y) {
        using Vector = Core::DenseVector<typename VectorType::ScalarType, VectorType::SizeAtCompile>;
        return spline(Vector::linspace(0, y.getLength() - 1, y.getLength()), y);
    }

    template<class VectorType1, class VectorType2>
    QSplineSeries& Plot::spline(const Core::RValueVector<VectorType1>& x, const Core::RValueVector<VectorType2>& y) {
        assert(x.getLength() == y.getLength());
        QSplineSeries* series = new QSplineSeries();
        for (size_t i = 0; i < x.getLength(); ++i)
            *series << QPointF(double(x.calc(i)), double(y.calc(i)));
        getChart()->addSeries(series);
        series->attachAxis(axisX);
        series->attachAxis(axisY);

        update();
        return *series;
    }

    template<class VectorType>
    QScatterSeries& Plot::scatter(const Core::RValueVector<VectorType>& y) {
        using Vector = Core::DenseVector<typename VectorType::ScalarType, VectorType::SizeAtCompile>;
        return scatter(Vector::linspace(0, y.getLength() - 1, y.getLength()), y);
    }

    template<class VectorType1, class VectorType2>
    QScatterSeries& Plot::scatter(const Core::RValueVector<VectorType1>& x, const Core::RValueVector<VectorType2>& y) {
        assert(x.getLength() == y.getLength());
        QScatterSeries* series = new QScatterSeries();
        for (size_t i = 0; i < x.getLength(); ++i)
            *series << QPointF(double(x.calc(i)), double(y.calc(i)));
        getChart()->addSeries(series);
        series->attachAxis(axisX);
        series->attachAxis(axisY);

        update();
        return *series;
    }

    template<class VectorType>
    QAreaSeries& Plot::hist(const Core::RValueVector<VectorType>& data, size_t binCount, bool density) {
        using ScalarType = typename VectorType::ScalarType;
        
        double binWidth, min;
        const size_t length = data.getLength();
        /* Get binWidth and min */ {
            ScalarType minimum = data.calc(0), maximum = data.calc(0);
            for (size_t i = 1; i < length; ++i) {
                ScalarType temp = data.calc(0);
                if (temp < minimum)
                    minimum = std::move(temp);
                else if (temp > maximum)
                    maximum = std::move(temp);
            }
            assert(maximum >= minimum);
            min = double(minimum);
            binWidth = double(maximum - minimum + (binCount - 1)) / binCount;
            if (binWidth == 0)
                binWidth = 1;
        }

        Core::Array<unsigned int> arr(binCount + 1, 0);
        const double binCountPerUnit = 1 / binWidth;
        for (size_t i = 0; i < length; ++i) {
            const size_t binIndex = size_t((double(data.calc(i)) - min) * binCountPerUnit);
            arr[binIndex]++;
        }

        QLineSeries* upper_series = new QLineSeries();
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
        QLineSeries* lower_series = new QLineSeries();
        *lower_series << QPointF(initial_x, 0) << QPointF(current_x, 0);

        QAreaSeries* series = new QAreaSeries(upper_series, lower_series);

        getChart()->addSeries(series);
        series->attachAxis(axisX);
        series->attachAxis(axisY);

        update();
        return *series;
    }

    template<class VectorType>
    QAreaSeries& Plot::area_boundary(const Core::RValueVector<VectorType>& x,
                                     const Core::RValueVector<VectorType>& lower,
                                     const Core::RValueVector<VectorType>& upper) {
        assert(x.getLength() == lower.getLength() && x.getLength() == upper.getLength());

        QLineSeries* upper_series = new QLineSeries();
        QLineSeries* lower_series = new QLineSeries();
        for (size_t i = 0; i < x.getLength(); ++i) {
            const double x_i = double(x.calc(i));
            *lower_series << QPointF(x_i, double(lower.calc(i)));
            *upper_series << QPointF(x_i, double(upper.calc(i)));
        }
        QAreaSeries* series = new QAreaSeries(upper_series, lower_series);
        getChart()->addSeries(series);
        series->attachAxis(axisX);
        series->attachAxis(axisY);
        update();
        return *series;
    }

    template<class VectorType>
    QAreaSeries& Plot::area_center(const Core::RValueVector<VectorType>& x,
                                   const Core::RValueVector<VectorType>& center,
                                   const Core::RValueVector<VectorType>& deviation) {
        const VectorType lower = center - deviation;
        const VectorType upper = center + deviation;
        return area_boundary(x, lower, upper);
    }

    template<class MatrixType>
    ContourSeries<MatrixType>& Plot::contour(const Core::LValueMatrix<MatrixType>& x,
                                             const Core::LValueMatrix<MatrixType>& y,
                                             const Core::LValueMatrix<MatrixType>& z,
                                             Core::Array<double> levels) {
        auto* series = new ContourSeries<MatrixType>(x, y, z, std::move(levels));
        series->attachTo(*getChart());
        series->attachAxis(axisX);
        series->attachAxis(axisY);
        update();
        return *series;
    }

    template<class VectorType>
    QBoxPlotSeries& Plot::boxWhisker(const Core::RValueVector<VectorType>& x, const Core::Array<VectorType>& data) {
        assert(x.getLength() == data.getLength());
        QBoxPlotSeries* series = new QBoxPlotSeries(QBoxPlotSeries::Numeric);
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

    template<class VectorType>
    QBoxPlotSeries& Plot::errorBar(
            const Core::RValueVector<VectorType>& mean,
            const Core::RValueVector<VectorType>& deviation) {
        using namespace Physica::Core;
        using ScalarType = float64;
        const auto x = VectorND<ScalarType>::linspace(0, mean.getLength() - 1, mean.getLength());
        return errorBar(x, mean, deviation);
    }

    template<class VectorType1, class VectorType2>
    QBoxPlotSeries& Plot::errorBar(
            const Core::RValueVector<VectorType1>& x,
            const Core::RValueVector<VectorType2>& mean,
            const Core::RValueVector<VectorType2>& deviation) {
        assert(x.getLength() == mean.getLength() && x.getLength() == deviation.getLength());
        QBoxPlotSeries* series = new QBoxPlotSeries(QBoxPlotSeries::Numeric);
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

    template<class VectorType>
    QBoxSet* Plot::setFromVector(const Core::LValueVector<VectorType>& v) {
        using BufferType = Core::DenseVector<typename VectorType::ScalarType, VectorType::SizeAtCompile>;
        BufferType buffer = v;
        std::sort(buffer.begin(), buffer.end());
        auto* result = new QBoxSet();
        const size_t count = v.getLength();
        result->setValue(QBoxSet::LowerExtreme, double(*buffer.begin()));
        result->setValue(QBoxSet::UpperExtreme, double(*buffer.rbegin()));
        result->setValue(QBoxSet::Median, findMedian(buffer, 0, count));
        result->setValue(QBoxSet::LowerQuartile, findMedian(buffer, 0, count / 2));
        result->setValue(QBoxSet::UpperQuartile, findMedian(buffer, count / 2 + (count % 2), count));
        return result;
    }

    template<class VectorType>
    double Plot::findMedian(const Core::LValueVector<VectorType>& sorted_v, size_t from, size_t to) {
        size_t count = to - from;
        if (count % 2) {
            return double(sorted_v[count / 2 + from]);
        }
        else {
            auto right = sorted_v[count / 2 + from];
            auto left = sorted_v[count / 2 - 1 + from];
            return double((right + left) * 0.5);
        }
    }

    inline void Plot::setAxisX(QValueAxis* axis) {
        auto* chart = Base::getChart();
        chart->removeAxis(axisX);
        chart->addAxis(axis, Qt::AlignBottom);
        axisX = axis;
    }

    inline void Plot::setAxisY(QValueAxis* axis) {
        auto* chart = Base::getChart();
        chart->removeAxis(axisY);
        chart->addAxis(axis, Qt::AlignLeft);
        axisY = axis;
    }

    inline void Plot::setAxisTop(QValueAxis* axis) {
        auto* chart = Base::getChart();
        chart->removeAxis(axisTop);
        chart->addAxis(axis, Qt::AlignTop);
        axisTop = axis;
    }

    inline void Plot::setAxisRight(QValueAxis* axis) {
        auto* chart = Base::getChart();
        chart->removeAxis(axisRight);
        chart->addAxis(axis, Qt::AlignRight);
        axisRight = axis;
    }

    inline void Plot::setMinX(double value) noexcept {
        axisX->setMin(value);
        axisTop->setMin(value);
    }

    inline void Plot::setMaxX(double value) noexcept {
        axisX->setMax(value);
        axisTop->setMax(value);
    }

    inline void Plot::setRangeX(double minX, double maxX) {
        setMinX(minX);
        setMaxX(maxX);
    }

    inline void Plot::setMinY(double value) noexcept {
        axisY->setMin(value);
        axisRight->setMin(value);
    }

    inline void Plot::setMaxY(double value) noexcept {
        axisY->setMax(value);
        axisRight->setMax(value);
    }

    inline void Plot::setRangeY(double minY, double maxY) {
        setMinY(minY);
        setMaxY(maxY);
    }

    inline void Plot::setDeltaX(double value) noexcept {
        axisX->setTickInterval(value);
        axisTop->setTickInterval(value);
    }

    inline void Plot::setDeltaY(double value) noexcept {
        axisY->setTickInterval(value);
        axisRight->setTickInterval(value);
    }
}
