/*
 * Copyright 2019-2021 WeiBo He.
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

#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QSplineSeries>
#include <QtCharts/QScatterSeries>
#include <QtCharts/QAreaSeries>
#include <Physica/Core/MultiPrecision/Scalar.h>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h>
#include "ContourSeries.h"

using Physica::Core::MultiScalar;

namespace Physica::Gui {
    class Plot : public QChartView {
    public:
        Plot(QWidget* parent = nullptr);
        Plot(MultiScalar (*func)(const MultiScalar&), const MultiScalar& begin
                , const MultiScalar& end, QWidget* parent = nullptr);
        /* Operations */
        template<class VectorType1, class VectorType2>
        QLineSeries& line(const Core::LValueVector<VectorType1>& x, const Core::LValueVector<VectorType2>& y);
        template<class VectorType>
        QSplineSeries& spline(const Core::LValueVector<VectorType>& y);
        template<class VectorType1, class VectorType2>
        QSplineSeries& spline(const Core::LValueVector<VectorType1>& x, const Core::LValueVector<VectorType2>& y);
        template<class VectorType>
        QScatterSeries& scatter(const Core::LValueVector<VectorType>& y);
        template<class VectorType1, class VectorType2>
        QScatterSeries& scatter(const Core::LValueVector<VectorType1>& x, const Core::LValueVector<VectorType2>& y);
        template<class VectorType>
        QAreaSeries& hist(const Core::LValueVector<VectorType>& data, size_t binCount, bool dencity = false);
        template<class MatrixType>
        ContourSeries<MatrixType>& contour(const Core::LValueMatrix<MatrixType>& x,
                                           const Core::LValueMatrix<MatrixType>& y,
                                           const Core::LValueMatrix<MatrixType>& z,
                                           Utils::Array<double> levels);
    };

    template<class VectorType1, class VectorType2>
    QLineSeries& Plot::line(const Core::LValueVector<VectorType1>& x, const Core::LValueVector<VectorType2>& y) {
        assert(x.getLength() == y.getLength());
        QLineSeries* series = new QLineSeries();
        for (size_t i = 0; i < x.getLength(); ++i)
            *series << QPointF(double(x[i]), double(y[i]));
        chart()->addSeries(series);

        update();
        return *series;
    }

    template<class VectorType>
    QSplineSeries& Plot::spline(const Core::LValueVector<VectorType>& y) {
        using Vector = Core::Vector<typename VectorType::ScalarType, VectorType::SizeAtCompile, VectorType::MaxSizeAtCompile>;
        return spline(Vector::linspace(0, y.getLength() - 1, y.getLength()), y);
    }

    template<class VectorType1, class VectorType2>
    QSplineSeries& Plot::spline(const Core::LValueVector<VectorType1>& x, const Core::LValueVector<VectorType2>& y) {
        assert(x.getLength() == y.getLength());
        QSplineSeries* series = new QSplineSeries();
        for (size_t i = 0; i < x.getLength(); ++i)
            *series << QPointF(double(x[i]), double(y[i]));
        chart()->addSeries(series);

        update();
        return *series;
    }

    template<class VectorType>
    QScatterSeries& Plot::scatter(const Core::LValueVector<VectorType>& y) {
        using Vector = Core::Vector<typename VectorType::ScalarType, VectorType::SizeAtCompile, VectorType::MaxSizeAtCompile>;
        return scatter(Vector::linspace(0, y.getLength() - 1, y.getLength()), y);
    }

    template<class VectorType1, class VectorType2>
    QScatterSeries& Plot::scatter(const Core::LValueVector<VectorType1>& x, const Core::LValueVector<VectorType2>& y) {
        assert(x.getLength() == y.getLength());
        QScatterSeries* series = new QScatterSeries();
        for (size_t i = 0; i < x.getLength(); ++i)
            *series << QPointF(double(x[i]), double(y[i]));
        chart()->addSeries(series);

        update();
        return *series;
    }

    template<class VectorType>
    QAreaSeries& Plot::hist(const Core::LValueVector<VectorType>& data, size_t binCount, bool dencity) {
        using ScalarType = typename VectorType::ScalarType;
        
        double binWidth, min;
        const size_t length = data.getLength();
        /* Get binWidth and min */ {
            ScalarType minimum = data[0], maximum = data[0];
            for (size_t i = 1; i < length; ++i) {
                ScalarType temp = data[i];
                if (temp < minimum)
                    minimum = std::move(temp);
                else if (temp > maximum)
                    maximum = std::move(temp);
            }
            assert(maximum >= minimum);
            min = double(minimum);
            binWidth = double(maximum - minimum) / binCount;
            if (binWidth == 0)
                binWidth = 1;
        }

        Utils::Array<unsigned int> arr(binCount + 1, 0);
        const double binCountPerUnit = 1 / binWidth;
        for (size_t i = 0; i < length; ++i) {
            const size_t binIndex = size_t((double(data[i]) - min) * binCountPerUnit);
            arr[binIndex]++;
        }

        QLineSeries* upper_series = new QLineSeries();
        double current_x = min;
        if (dencity) {
            const double dencity_factor = 1 / (binWidth * length);
            for (size_t i = 0; i < binCount; ++i) {
                const double y = arr[i] * dencity_factor;
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
        *lower_series << QPointF(min, 0) << QPointF(current_x, 0);

        QAreaSeries* series = new QAreaSeries(upper_series, lower_series);

        chart()->addSeries(series);

        update();
        return *series;
    }

    template<class MatrixType>
    ContourSeries<MatrixType>& Plot::contour(const Core::LValueMatrix<MatrixType>& x,
                                             const Core::LValueMatrix<MatrixType>& y,
                                             const Core::LValueMatrix<MatrixType>& z,
                                             Utils::Array<double> levels) {
        auto* series = new ContourSeries<MatrixType>(x, y, z, std::move(levels));
        series->attachTo(*chart());
        update();
        return *series;
    }
}
