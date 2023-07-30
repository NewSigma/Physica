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
#include "Physica/Core/Physics/Phonon/FinitePhonon.h"

namespace Physica::Gui {
    template<class ScalarType, class PosScalarType>
    class PhononPlot : public Plot {
        using Vector3D = Core::Vector<ScalarType, 3>;
        using VectorType = Core::Vector<ScalarType>;
        using PhononType = Core::FinitePhonon<ScalarType, PosScalarType>;
        using MatrixGrid = typename PhononType::MatrixGrid;

        ScalarType currentX;
    public:
        PhononPlot();
        /* Operations */
        void plotPathScatter(
            const PhononType& phonon,
            const MatrixGrid& forceConstants,
            Vector3D from,
            Vector3D to,
            size_t numPoint,
            const char* label);
        void plotPathLine(
            const PhononType& phonon,
            const MatrixGrid& forceConstants,
            Vector3D from,
            Vector3D to,
            size_t numPoint,
            const char* label);
    };

    template<class ScalarType, class PosScalarType>
    PhononPlot<ScalarType, PosScalarType>::PhononPlot() : Plot(0, 0, 0, 0, 1, 100), currentX(0) {
        const QFont font = Plot::getAxisX()->labelsFont();
        QCategoryAxis* axisX = new QCategoryAxis();
        axisX->append("Γ", 0);
        axisX->setLabelsPosition(QCategoryAxis::AxisLabelsPositionOnValue);
        axisX->setLinePenColor(Qt::black);
        axisX->setLabelsFont(font);
        Plot::setAxisX(axisX);

        Plot::getAxisY()->setLabelFormat("%d");
        Plot::getAxisY()->setTitleText("Frequency/THz");
    }

    template<class ScalarType, class PosScalarType>
    void PhononPlot<ScalarType, PosScalarType>::plotPathScatter(
            const PhononType& phonon,
            const MatrixGrid& forceConstants,
            Vector3D from,
            Vector3D to,
            size_t numPoint,
            const char* label) {
        using namespace Physica::Core;
        const size_t numBranch = phonon.getUnitCellDOF();
        const auto factors = VectorType::linspace(0, 1, numPoint);
        const ScalarType deltaX = (to - from).norm();
        const VectorType x = currentX + factors * deltaX;

        DenseMatrix<ScalarType> buffer(numBranch, numPoint);
        ScalarType minFreq = Plot::getMinY();
        ScalarType maxFreq = Plot::getMaxY();
        for (size_t i = 0; i < numPoint; ++i) {
            const ScalarType factor = factors[i];
            const Vector3D qPoint = from * (ScalarType(1) - factor) + to * factor;
            auto fcMatrix = phonon.interpolatePoint(qPoint, forceConstants);
            phonon.toDynamicMatrix(fcMatrix);
            const auto eigen = PhononType::diagonalize(fcMatrix);
            auto freq = phonon.makeFreq(eigen);
            freq *= ScalarType(1E-12 / PhyConst<AU>::timeToSecond(1));
            buffer.col(i) = freq;
            minFreq = std::min(minFreq, freq.min());
            maxFreq = std::max(maxFreq, freq.max());
        }

        for (size_t i = 0; i < numBranch; ++i) {
            auto& scatter = Plot::scatter(x, buffer.row(i));
            auto pen = scatter.pen();
            pen.setColor(Qt::black);
            scatter.setPen(pen);
            scatter.setMarkerSize(2);
            scatter.setColor(Qt::black);
        }
        Plot::setMinY(double(minFreq));
        Plot::setMaxY(double(maxFreq));
        currentX += deltaX;
        Plot::setMaxX(double(currentX));
        static_cast<QCategoryAxis*>(Plot::getAxisX())->append(label, double(currentX));
    }
    /**
     * Estimating path connection algorithm from phonopy[1]
     * Reference:
     * [1] https://github.com/phonopy/phonopy
     */
    template<class ScalarType, class PosScalarType>
    void PhononPlot<ScalarType, PosScalarType>::plotPathLine(
            const PhononType& phonon,
            const MatrixGrid& forceConstants,
            Vector3D from,
            Vector3D to,
            size_t numPoint,
            const char* label) {
        using namespace Physica::Core;
        const size_t numBranch = phonon.getUnitCellDOF();
        const auto factors = VectorType::linspace(0, 1, numPoint);
        const ScalarType deltaX = (to - from).norm();
        const VectorType x = currentX + factors * deltaX;

        DenseMatrix<ScalarType> buffer(numBranch, numPoint);
        DenseMatrix<ScalarType> lastEigenvectors;
        ScalarType minFreq = Plot::getMinY();
        ScalarType maxFreq = Plot::getMaxY();
        for (size_t i = 0; i < numPoint; ++i) {
            const ScalarType factor = factors[i];
            const Vector3D qPoint = from * (ScalarType(1) - factor) + to * factor;
            auto fcMatrix = phonon.interpolatePoint(qPoint, forceConstants);
            phonon.toDynamicMatrix(fcMatrix);
            const auto eigen = PhononType::diagonalize(fcMatrix);
            auto freq = phonon.makeFreq(eigen);
            freq *= ScalarType(1E-12 / PhyConst<AU>::timeToSecond(1));
            minFreq = std::min(minFreq, freq.min());
            maxFreq = std::max(maxFreq, freq.max());

            const bool shouldInit = i == 0;
            auto eigenvectors = phonon.makeEigenVectors(eigen);
            auto bufferCol = buffer.col(i);
            if (shouldInit) {
                bufferCol = freq;
                lastEigenvectors.swap(eigenvectors);
            }
            else {
                Utils::Array<size_t> occupiedIndex{};
                occupiedIndex.reserve(numBranch);
                for (size_t branch = 0; branch < numBranch; ++branch) {
                    const VectorType dots = abs(eigenvectors.transpose() * lastEigenvectors.col(branch));
                    size_t index = numBranch;
                    ScalarType maxDot = 0;
                    for (size_t j = 0; j < dots.getLength(); ++j) {
                        const ScalarType dot = abs(dots[j]);
                        const bool isNotOccupied = std::find(occupiedIndex.cbegin(), occupiedIndex.cend(), j) == occupiedIndex.cend();
                        if (isNotOccupied && dot > maxDot) {
                            index = j;
                            maxDot = dot;
                        }
                    }
                    assert(index != numBranch && "[Error]: Unexpected result. Something wrong with the dynamic matrix?");
                    occupiedIndex.append(index);
                    bufferCol[branch] = freq[index];
                    lastEigenvectors.asArray()[branch].swap(eigenvectors.asArray()[index]);
                }
            }
        }

        for (size_t i = 0; i < numBranch; ++i) {
            auto& line = Plot::line(x, buffer.row(i));
            line.setColor(Qt::black);
        }
        Plot::setMinY(double(minFreq));
        Plot::setMaxY(double(maxFreq));
        currentX += deltaX;
        Plot::setMaxX(double(currentX));
        static_cast<QCategoryAxis*>(Plot::getAxisX())->append(label, double(currentX));
    }
}
