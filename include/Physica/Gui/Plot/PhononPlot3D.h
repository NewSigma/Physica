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

#include "Physica/Core/Physics/Phonon/FrozenPhonon.h"
#include "Plot3D.h"

namespace Physica::Gui {
    template<Scalar T>
    class PhononPlot3D : public Plot3D {
        using Base = Plot3D;
        using Vector3D = Core::Vector3D<T>;
        using VectorType = Core::VectorND<T>;
        using MatrixType = Core::DenseMatrix<T>;
        using MeshType = std::pair<MatrixType, MatrixType>;
        using BandArray = Core::Array<MatrixType>;
        using PhononType = Core::FrozenPhonon<T>;
        using KSpaceFCGrid = PhononType::KSpaceFCGrid;
    public:
        PhononPlot3D() = default;
        PhononPlot3D(const PhononPlot3D&) = default;
        PhononPlot3D(PhononPlot3D&&) noexcept = default;
        ~PhononPlot3D() = default;
        /* Operators */
        PhononPlot3D& operator=(const PhononPlot3D&) = default;
        PhononPlot3D& operator=(PhononPlot3D&&) noexcept = default;
        /* Operations */
        QList<QSurface3DSeries*> plotBand(const PhononType& ph, const KSpaceFCGrid& forceConstants, size_t numPointX, size_t numPointY);
    private:
        BandArray calcBands(const PhononType& ph, const KSpaceFCGrid& forceConstants, const MeshType& mesh);
    };

    template<Scalar T>
    QList<QSurface3DSeries*> PhononPlot3D<T>::plotBand(
            const PhononType& ph, const KSpaceFCGrid& forceConstants, size_t numPointX, size_t numPointY) {
        const auto x = VectorType::linspace(0, 1, numPointX);
        const auto y = VectorType::linspace(0, 1, numPointY);
        const MeshType mesh = MatrixType::meshgrid(x, y);
        const auto bands = calcBands(ph, forceConstants, mesh);
        const QLinearGradient gr = Base::makeDefaultGrad();
        for (size_t i = 0; i < bands.getLength(); ++i) {
            auto& surface = Base::surf(mesh.first, mesh.second, bands[i]);
            surface.setBaseGradient(gr);
            surface.setColorStyle(Q3DTheme::ColorStyleRangeGradient);
        }

        auto& axisX = *Base::getAxisX();
        axisX.setRange(0, 1);
        axisX.setLabelFormat("%.1f");
        axisX.setTitle("X");
        axisX.setTitleVisible(true);

        auto& axisY = *Base::getAxisY();
        axisY.setRange(0, 1);
        axisY.setLabelFormat("%.1f");
        axisY.setTitle("Y");
        axisY.setTitleVisible(true);

        auto& axisZ = *Base::getAxisZ();
        axisZ.setLabelFormat("%d");
        axisZ.setTitle("Frequency/THz");
        axisZ.setTitleVisible(true);
        return Base::surface->seriesList();
    }

    template<Scalar T>
    PhononPlot3D<T>::BandArray PhononPlot3D<T>::calcBands(
            const PhononType& ph, const KSpaceFCGrid& forceConstants, const MeshType& mesh) {
        using namespace Physica::Core;
        const MatrixType& meshX = mesh.first;
        const MatrixType& meshY = mesh.second;
        BandArray result(ph.getNumBand(), meshX.getRow(), meshX.getCol());
        T minFreq = Base::getMinZ();
        T maxFreq = Base::getMaxZ();
        for (size_t major = 0; major < meshX.getMaxMajor(); ++major) {
            for (size_t minor = 0; minor < meshX.getMaxMinor(); ++minor) {
                const Vector3D qPoint{meshX.calcFromMajorMinor(major, minor), meshY.calcFromMajorMinor(major, minor), 0};
                auto fcMatrix = ph.interpolatePoint(qPoint, forceConstants);
                ph.toDynamicMatrix(fcMatrix);
                auto eigen = PhononType::diagonalize(fcMatrix);
                auto freq = ph.makeFreq(eigen);
                freq *= T(PhyConst<AU>::freqToTHz(1));

                minFreq = std::min(minFreq, freq.min());
                maxFreq = std::max(maxFreq, freq.max());
                for (size_t i = 0; i < result.getLength(); ++i)
                    result[i].refFromMajorMinor(major, minor) = freq[i];
            }
        }
        Base::setMinZ(float(minFreq));
        Base::setMaxZ(float(maxFreq));
        return result;
    }
}
