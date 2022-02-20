/*
 * Copyright 2021-2022 WeiBo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "KPoint.h"
#include "Grid3D.h"

namespace Physica::Core {
    template<class ScalarType, bool isSpinPolarized> class BandGrid;

    template<class ScalarType, bool isSpinPolarized>
    class BandGrid {
        using KPoints = Grid3D<KPoint<ScalarType, isSpinPolarized>, false>;
        using LatticeMatrix = typename KPoints::LatticeMatrix;

        KPoints kPoints;
        size_t electronCount;
    public:
        template<class MatrixType>
        BandGrid(ScalarType cutEnergy, const LValueMatrix<MatrixType>& repLatt, size_t kPointX, size_t kPointY, size_t kPointZ, size_t electronCount_);
        BandGrid(const BandGrid&) = default;
        BandGrid(BandGrid&&) noexcept = default;
        ~BandGrid() = default;
        /* Operators */
        BandGrid& operator=(BandGrid band) noexcept;
        /* Getters */
        [[nodiscard]] KPoints& getKPoints() { return kPoints; }
        [[nodiscard]] const KPoints& getKPoints() const noexcept { return kPoints; }
        [[nodiscard]] ScalarType getTotalEnergy() const noexcept;
        template<class VectorType>
        [[nodiscard]] VectorType getDensityOfStates(const LValueVector<VectorType>& atEnergy) const;
        /* Helpers */
        void swap(BandGrid& band) noexcept;
    private:
        [[nodiscard]] Vector<ScalarType, 3> gradEnergy(size_t kPointId) const;
    };

    template<class ScalarType, bool isSpinPolarized>
    template<class MatrixType>
    BandGrid<ScalarType, isSpinPolarized>::BandGrid(ScalarType cutEnergy,
                                                    const LValueMatrix<MatrixType>& repLatt,
                                                    size_t kPointX,
                                                    size_t kPointY,
                                                    size_t kPointZ,
                                                    size_t electronCount_)
            : kPoints(repLatt, kPointX, kPointY, kPointZ)
            , electronCount(electronCount_) {
        assert(kPoints.getSize() != 0);
        size_t kPointID = 0;
        const size_t plainWaveCount = Grid3D<double, true>::sizeFromCutEnergy(cutEnergy, repLatt); //TODO: signed/unsigned character can be moved to father class
        
        const ScalarType kPointWeight = reciprocal(ScalarType(kPoints.getSize()));
        const ScalarType stepX = reciprocal(ScalarType(kPointX));
        const ScalarType stepY = reciprocal(ScalarType(kPointY));
        const ScalarType stepZ = reciprocal(ScalarType(kPointZ));

        Vector<ScalarType, 3> k{};
        ScalarType& kx = k[0];
        ScalarType& ky = k[1];
        ScalarType& kz = k[2];

        kx = (ScalarType(1) - ScalarType(kPointX)) / ScalarType(2 * kPointX);
        for (size_t x = 1; x <= kPointX; ++x) {
            ky = (ScalarType(1) - ScalarType(kPointY)) / ScalarType(2 * kPointY);
            for (size_t y = 1; y <= kPointY; ++y) {
                kz = (ScalarType(1) - ScalarType(kPointZ)) / ScalarType(2 * kPointZ);
                for (size_t z = 1; z <= kPointZ; ++z) {
                    kPoints[kPointID] = KPoint<ScalarType, isSpinPolarized>(k, plainWaveCount, kPointWeight);
                    kz += stepZ;
                    ++kPointID;
                }
                ky += stepY;
            }
            kx += stepX;
        }
    }

    template<class ScalarType, bool isSpinPolarized>
    ScalarType BandGrid<ScalarType, isSpinPolarized>::getTotalEnergy() const noexcept {
        ScalarType energy = ScalarType::Zero();
        for (auto ite = kPoints.cbegin(); ite != kPoints.cend(); ++ite)
            energy += (*ite).getTotalEnergy();
        energy *= reciprocal(ScalarType(kPoints.getSize()));
        return energy;
    }
    /**
     * Reference:
     * [1] Bross H. On the Efficiency of Different Schemes for the Evaluation of the Density of States and Related Properties in Solids[J]. Physica Status Solidi, 2010, 179(2):429-439.
     */
    template<class ScalarType, bool isSpinPolarized>
    template<class VectorType>
    VectorType BandGrid<ScalarType, isSpinPolarized>::getDensityOfStates(const LValueVector<VectorType>& atEnergy) const {
        auto dos = VectorType(atEnergy.getLength());
        for (size_t i = 0; i < atEnergy.getLength(); ++i) {
            const ScalarType energy = atEnergy[i];
            ScalarType density = ScalarType::Zero();
            for (size_t kPointId = 0; kPointId < kPoints.getSize(); ++i) {
                const ScalarType energy0 = kPoints[kPointId].getTotalEnergy();
                const auto gradE = gradEnergy(kPointId);
                const ScalarType normalizer = gradE[0] * gradE[1] * gradE[2] * ScalarType(0.5);

                ScalarType deltaDensity = ScalarType::Zero();
                for (int sigma1; sigma1 < 2; ++sigma1) {
                    for (int sigma2; sigma2 < 2; ++sigma2) {
                        for (int sigma3; sigma3 < 2; ++sigma3) {
                            const int sum = sigma1 + sigma2 + sigma3;
                            ScalarType temp = energy - energy0;
                            temp -= ScalarType(sigma1 == 0 ? 1 : -1) * gradE(0);
                            temp -= ScalarType(sigma2 == 0 ? 1 : -2) * gradE(1);
                            temp -= ScalarType(sigma3 == 0 ? 1 : -3) * gradE(2);
                            if (temp.isPositive()) {
                                temp = square(temp);
                                deltaDensity += (sum % 2 == 0 ? temp : -temp);
                            }
                        }
                    }
                }
                density += deltaDensity * normalizer;
            }
            dos[i] = density;
        }
        return dos;
    }

    template<class ScalarType, bool isSpinPolarized>
    BandGrid<ScalarType, isSpinPolarized>& BandGrid<ScalarType, isSpinPolarized>::operator=(BandGrid band) noexcept {
        swap(band);
        return *this;
    }

    template<class ScalarType, bool isSpinPolarized>
    void BandGrid<ScalarType, isSpinPolarized>::swap(BandGrid& band) noexcept {
        swap(kPoints, band.kPoints);
        swap(electronCount, band.electronCount);
    }

    template<class ScalarType, bool isSpinPolarized>
    inline void swap(Physica::Core::BandGrid<ScalarType, isSpinPolarized>& band1,
                     Physica::Core::BandGrid<ScalarType, isSpinPolarized>& band2) noexcept {
        band1.swap(band2);
    }
    /**
     * \returns Gradient of energy in t-coordinate defined in [1]
     * 
     * Reference:
     * [1] Bross H. On the Efficiency of Different Schemes for the Evaluation of the Density of States and Related Properties in Solids[J]. Physica Status Solidi, 2010, 179(2):429-439.
     */
    template<class ScalarType, bool isSpinPolarized>
    Vector<ScalarType, 3> BandGrid<ScalarType, isSpinPolarized>::gradEnergy(size_t kPointId) const {
        auto dimAdd = [](size_t dim, size_t dim_all) { return dim == dim_all - 1 ? 0 : dim + 1; };
        auto dimSub = [](size_t dim, size_t dim_all) { return dim == 0 ? dim_all - 1 : dim - 1; };
        const auto[x, y, z] = kPoints.indexToDim(kPointId);
        const size_t dimX = kPoints.getDimX();
        const size_t dimY = kPoints.getDimY();
        const size_t dimZ = kPoints.getDimZ();
        const ScalarType factor = ScalarType(0.25);
        const ScalarType gradX = (kPoints(dimAdd(x, dimX), y, z) -  kPoints(dimSub(x, dimX), y, z)) * factor;
        const ScalarType gradY = (kPoints(x, dimAdd(y, dimY), z) -  kPoints(x, dimSub(y, dimY), z)) * factor;
        const ScalarType gradZ = (kPoints(x, y, dimAdd(z, dimZ)) -  kPoints(x, y, dimSub(z, dimZ))) * factor;
        return {gradX, gradY, gradZ};
    }
}
