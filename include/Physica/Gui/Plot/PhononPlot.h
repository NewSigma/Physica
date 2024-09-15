/*
 * Copyright 2023-2024 Weibo He.
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
        using MatrixType = Core::DenseMatrix<ScalarType>;
        using ComplexMatrix = Core::DenseMatrix<ComplexType>;
        using PhononType = Core::FrozenPhonon<ScalarType>;
        using KSpaceFCGrid = typename PhononType::KSpaceFCGrid;
        using EigenSolverType = typename PhononType::EigenSolverType;
        using Index3D = typename Core::GridBase::Index3D;
        using ColorArray = Core::Array<QColor>;

        enum BandConnectMethod {
            Direct,
            Predict
        };
    private:
        ScalarType currentX;
        ColorArray colors;
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
        void plotBandLine(
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
        /* Setters */
        void resetCurrentX() { currentX = ScalarType(0); }
    private:
        EigenSolverType diagonalize(
            const PhononType& phonon,
            const KSpaceFCGrid& forceConstants,
            ScalarType factor,
            Vector3D from,
            Vector3D to);
        VectorType makeFreq(const EigenSolverType& eigen);
        MatrixType calcPaths(
            const PhononType& phonon,
            const KSpaceFCGrid& forceConstants,
            const VectorType& factors,
            Vector3D from,
            Vector3D to,
            BandConnectMethod method);
        void predictConnect(
            size_t point,
            VectorType& freq,
            ComplexMatrix* eigenvector,
            const VectorType& lastFreq,
            VectorType& lastDeltaFreq);
        void sortFreq(VectorType& freq, ComplexMatrix* eigenvector, const VectorType& predictFreq) const;
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
        currentX += deltaX;
        Plot::setMaxX(double(currentX));
        getAxisX()->append(label, double(currentX));
    }

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
        const auto paths = calcPaths(phonon, forceConstants, factors, from, to, method);
        for (size_t i = 0; i < phonon.getNumBand(); ++i) {
            const auto freq = paths.row(i);
            auto& line = Plot::line(x, freq);
            line.setColor(Qt::black);
        }
        currentX += deltaX;
        Plot::setMaxX(double(currentX));
        getAxisX()->append(label, double(currentX));
    }

    template<class ScalarType>
    void PhononPlot<ScalarType>::plotBandLine(
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
        const size_t numBand = phonon.getNumBand();
        colors.resize(numBand);

        VectorType freq{}, lastFreq{}, lastDeltaFreq{}, bufferX{}, bufferY{};
        ComplexMatrix eigenvector{}, lastEigenvector{};
        ColorArray colorBuffer = colors;
        for (size_t p = 0; p < numPoint; ++p) {
            lastFreq.swap(freq);
            lastEigenvector.swap(eigenvector);

            const auto eigen = diagonalize(phonon, forceConstants, factors[p], from, to);
            freq = makeFreq(eigen);
            eigenvector = eigen.getEigenvectors();
            if (method == BandConnectMethod::Predict)
                predictConnect(p, freq, &eigenvector, lastFreq, lastDeltaFreq);
            const bool isLastFreqReady = p > 0;
            if (!isLastFreqReady)
                continue;

            for (size_t band = 0; band < numBand; ++band) {
                bufferX.append(x[p - 1]);
                bufferY.append(lastFreq[band]);
                for (size_t otherBand = 0; otherBand < numBand; ++otherBand) {
                    if (band == otherBand)
                        continue;
                    const ScalarType deltaFreq1 = lastFreq[band] - lastFreq[otherBand];
                    const ScalarType deltaFreq2 = freq[band] - freq[otherBand];
                    const bool noCross = (deltaFreq1 * deltaFreq2).isPositive();
                    if (noCross)
                        continue;

                    const ScalarType dot11 = (lastEigenvector.col(band).conjugate() * eigenvector.col(band)).norm();
                    const ScalarType dot12 = (lastEigenvector.col(band).conjugate() * eigenvector.col(otherBand)).norm();
                    const ScalarType dot21 = (lastEigenvector.col(otherBand).conjugate() * eigenvector.col(band)).norm();
                    const ScalarType dot22 = (lastEigenvector.col(otherBand).conjugate() * eigenvector.col(otherBand)).norm();
                    const bool isReversion = (dot12 > dot11) && (dot21 > dot22);
                    if (isReversion) {
                        const ScalarType deltaFreq3 = lastFreq[band] - freq[band];
                        const ScalarType deltaFreq4 = lastFreq[otherBand] - freq[otherBand];
                        const ScalarType factor = reciprocal(deltaFreq3 - deltaFreq4);
                        const ScalarType crossX = (-deltaFreq2 * x[p - 1] + deltaFreq1 * x[p]) * factor;
                        const ScalarType crossY = (lastFreq[band] * freq[otherBand] - lastFreq[otherBand] * freq[band]) * factor;
                        bufferX.append(crossX);
                        bufferY.append(crossY);
                        bufferX.append(x[p]);
                        bufferY.append(freq[otherBand]);
                        if (band < otherBand)
                            std::swap(colorBuffer[band], colorBuffer[otherBand]);
                        goto done;
                    }
                }
                bufferX.append(x[p]);
                bufferY.append(freq[band]);

            done:
                auto& line = Plot::line(bufferX, bufferY);
                if (currentX.isZero() && p == 1)
                    colorBuffer[band] = line.color();
                else
                    line.setColor(colors[band]);
                bufferX.clear();
                bufferY.clear();
            }
            colors = colorBuffer;
        }
        currentX += deltaX;
        Plot::setMaxX(double(currentX));
        getAxisX()->append(label, double(currentX));
    }

    template<class ScalarType>
    typename PhononPlot<ScalarType>::EigenSolverType PhononPlot<ScalarType>::diagonalize(
            const PhononType& phonon,
            const KSpaceFCGrid& forceConstants,
            ScalarType factor,
            Vector3D from,
            Vector3D to) {
        const Vector3D qPoint = from * (ScalarType(1) - factor) + to * factor;
        auto fcMatrix = phonon.interpolatePoint(qPoint, forceConstants);
        phonon.toDynamicMatrix(fcMatrix);
        return PhononType::diagonalize(fcMatrix);
    }

    template<class ScalarType>
    typename PhononPlot<ScalarType>::VectorType
    PhononPlot<ScalarType>::makeFreq(const EigenSolverType& eigen) {
        using namespace Physica::Core;
        auto freq = PhononType::makeFreq(eigen);
        freq *= ScalarType(Core::PhyConst<AU>::freqToTHz(1));
        Plot::setMinY(std::min(Plot::getMinY(), double(freq.min())));
        Plot::setMaxY(std::max(Plot::getMaxY(), double(freq.max())));
        return freq;
    }

    template<class ScalarType>
    typename PhononPlot<ScalarType>::MatrixType PhononPlot<ScalarType>::calcPaths(
            const PhononType& phonon,
            const KSpaceFCGrid& forceConstants,
            const VectorType& factors,
            Vector3D from,
            Vector3D to,
            BandConnectMethod method) {
        const size_t numPoint = factors.getLength();
        MatrixType result(phonon.getNumBand(), numPoint);
        VectorType lastDeltaFreq;
        for (size_t i = 0; i < numPoint; ++i) {
            const auto eigen = diagonalize(phonon, forceConstants, factors[i], from, to);
            auto freq = makeFreq(eigen);
            if (method == BandConnectMethod::Predict)
                predictConnect(i, freq, nullptr, result.asArray()[i == 0 ? 0 : (i - 1)], lastDeltaFreq);
            auto col = result.col(i);
            col = freq;
        }
        return result;
    }

    template<class ScalarType>
    void PhononPlot<ScalarType>::predictConnect(
            size_t point,
            VectorType& freq,
            ComplexMatrix* eigenvector,
            const VectorType& lastFreq,
            VectorType& lastDeltaFreq) {
        const bool isLastFreqReady = point != 0;
        if (isLastFreqReady) {
            const bool isLastDeltaFreqReady = point > 1;
            if (isLastDeltaFreqReady) {
                const VectorType predictFreq = lastFreq + lastDeltaFreq;
                sortFreq(freq, eigenvector, predictFreq);
            }
            lastDeltaFreq = freq - lastFreq;
        }
    }

    template<class ScalarType>
    void PhononPlot<ScalarType>::sortFreq(VectorType& freq, ComplexMatrix* eigenvector, const VectorType& predictFreq) const {
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
            if (index != band) {
                freq[index].swap(freq[band]);
                if (eigenvector != nullptr) {
                    auto& arr = (*eigenvector).asArray();
                    arr[index].swap(arr[band]);
                }
            }
        }
    }
}
