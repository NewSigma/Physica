/*
 * Copyright 2023 WeiBo He.
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

#include <QtCharts/QCategoryAxis>
#include "Plot.h"
#include "Physica/Core/Physics/Phonon/FrozenPhonon.h"
#include "Physica/Core/Physics/SolidState/BrillouInterpolate.h"

namespace Physica::Gui {
    template<class ScalarType>
    class PhononPlot : public Plot {
    public:
        using ComplexType = Core::ComplexScalar<ScalarType>;
        using Vector3D = Core::Vector<ScalarType, 3>;
        using VectorType = Core::Vector<ScalarType>;
        using PhononType = Core::FrozenPhonon<ScalarType>;
        using KSpaceFCGrid = typename PhononType::KSpaceFCGrid;
        using Index3D = typename Core::GridBase::Index3D;

        enum BandConnectMethod {
            Direct,
            Predict
        };
    private:
        ScalarType currentX;
    public:
        PhononPlot();
        /* Operations */
        void plotPathScatter(
            const PhononType& phonon,
            const KSpaceFCGrid& forceConstants,
            Vector3D from,
            Vector3D to,
            size_t numPoint,
            const char* label);
        void plotPathLine(
            const PhononType& phonon,
            const KSpaceFCGrid& forceConstants,
            Vector3D from,
            Vector3D to,
            size_t numPoint,
            const char* label,
            BandConnectMethod method = BandConnectMethod::Predict);
        /* Getters */
        [[nodiscard]] QCategoryAxis* getAxisX() const noexcept { return reinterpret_cast<QCategoryAxis*>(Plot::getAxisX()); }
        [[nodiscard]] ScalarType getMaxFreqInAU() const noexcept { return getMaxY() / ScalarType(Core::PhyConst<Core::AU>::freqToTHz(1)); }
        [[nodiscard]] ScalarType getMinFreqInAU() const noexcept { return getMinY() / ScalarType(Core::PhyConst<Core::AU>::freqToTHz(1)); }
    private:
        Core::DenseMatrix<ScalarType> calcPaths(
            const PhononType& phonon,
            const KSpaceFCGrid& forceConstants,
            const VectorType& factors,
            Vector3D from,
            Vector3D to,
            BandConnectMethod method);
        void sortFreq(VectorType& freq, const VectorType& predictFreq) const;
    };

    template<class ScalarType>
    PhononPlot<ScalarType>::PhononPlot() : Plot(0, 0, 0, 0, 1, 100), currentX(0) {
        const QFont font = Plot::getAxisX()->labelsFont();
        QCategoryAxis* axisX = new QCategoryAxis();
        axisX->append("Γ", 0);
        axisX->setLabelsPosition(QCategoryAxis::AxisLabelsPositionOnValue);
        axisX->setLinePenColor(Qt::black);
        axisX->setLabelsFont(font);
        axisX->setTitleFont(font);
        Plot::setAxisX(axisX);

        Plot::getAxisY()->setLabelFormat("%d");
        Plot::getAxisY()->setTitleText("Frequency/THz");
    }

    template<class ScalarType>
    void PhononPlot<ScalarType>::plotPathScatter(
            const PhononType& phonon,
            const KSpaceFCGrid& forceConstants,
            Vector3D from,
            Vector3D to,
            size_t numPoint,
            const char* label) {
        const auto factors = VectorType::linspace(0, 1, numPoint);
        const ScalarType deltaX = (to - from).norm();
        const VectorType x = currentX + factors * deltaX;
        currentX += deltaX;
        Plot::setMaxX(double(currentX));

        const auto paths = calcPaths(phonon, forceConstants, factors, from, to, BandConnectMethod::Direct);
        for (size_t i = 0; i < phonon.getNumBand(); ++i) {
            const auto freq = paths.row(i);
            auto& scatter = Plot::scatter(x, freq);
            auto pen = scatter.pen();
            pen.setColor(Qt::black);
            scatter.setPen(pen);
            scatter.setMarkerSize(2);
            scatter.setColor(Qt::black);
        }
        getAxisX()->append(label, double(currentX));
    }
    /**
     * Estimating path connection algorithm from phonopy[1]
     * Reference:
     * [1] https://github.com/phonopy/phonopy
     */
    template<class ScalarType>
    void PhononPlot<ScalarType>::plotPathLine(
            const PhononType& phonon,
            const KSpaceFCGrid& forceConstants,
            Vector3D from,
            Vector3D to,
            size_t numPoint,
            const char* label,
            BandConnectMethod method) {
        const auto factors = VectorType::linspace(0, 1, numPoint);
        const ScalarType deltaX = (to - from).norm();
        const VectorType x = currentX + factors * deltaX;
        currentX += deltaX;
        Plot::setMaxX(double(currentX));

        const auto paths = calcPaths(phonon, forceConstants, factors, from, to, method);
        for (size_t i = 0; i < phonon.getNumBand(); ++i) {
            const auto freq = paths.row(i);
            auto& line = Plot::line(x, freq);
            line.setColor(Qt::black);
        }
        getAxisX()->append(label, double(currentX));
    }

    template<class ScalarType>
    Core::DenseMatrix<ScalarType> PhononPlot<ScalarType>::calcPaths(
            const PhononType& phonon,
            const KSpaceFCGrid& forceConstants,
            const VectorType& factors,
            Vector3D from,
            Vector3D to,
            BandConnectMethod method) {
        using namespace Physica::Core;
        const size_t numPoint = factors.getLength();
        DenseMatrix<ScalarType> result(phonon.getNumBand(), numPoint);
        ScalarType minFreq = Plot::getMinY();
        ScalarType maxFreq = Plot::getMaxY();

        VectorType lastFreq, lastDeltaFreq;
        for (size_t i = 0; i < numPoint; ++i) {
            const ScalarType factor = factors[i];
            const Vector3D qPoint = from * (ScalarType(1) - factor) + to * factor;
            auto fcMatrix = phonon.interpolatePoint(qPoint, forceConstants);
            phonon.toDynamicMatrix(fcMatrix);
            auto eigen = PhononType::diagonalize(fcMatrix);
            auto freq = phonon.makeFreq(eigen);
            freq *= ScalarType(PhyConst<AU>::freqToTHz(1));
            if (method == BandConnectMethod::Predict) {
                const bool isLastFreqReady = i != 0;
                if (isLastFreqReady != 0) {
                    const VectorType predictFreq = lastFreq + lastDeltaFreq;
                    const bool isLastDeltaFreqReady = i > 1;
                    if (isLastDeltaFreqReady)
                        sortFreq(freq, predictFreq);
                    lastDeltaFreq = freq - lastFreq;
                }
                lastFreq = freq;
            }
            auto col = result.col(i);
            col = freq;
            minFreq = std::min(minFreq, freq.min());
            maxFreq = std::max(maxFreq, freq.max());
        }
        Plot::setMinY(double(minFreq));
        Plot::setMaxY(double(maxFreq));
        return result;
    }

    template<class ScalarType>
    void PhononPlot<ScalarType>::sortFreq(VectorType& freq, const VectorType& predictFreq) const {
        const size_t numBand = predictFreq.getLength();
        for (size_t band = 0; band < numBand - 1; ++band) {
            ScalarType minDelta = std::numeric_limits<ScalarType>::max();
            size_t index = band;
            for (size_t otherBand = band; otherBand < numBand; ++otherBand) {
                const ScalarType delta = abs(freq[otherBand] - predictFreq[band]);
                if (delta < minDelta) {
                    minDelta = delta;
                    index = otherBand;
                }
            }
            if (index != band)
                freq[index].swap(freq[band]);
        }
    }
}
