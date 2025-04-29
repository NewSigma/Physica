/*
 * Copyright 2021-2022 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/Core/Utils/Container/ArrayND.h"
#include "KPoint.h"

namespace Physica {
    template<Scalar T, bool isSpinPolarized> class BandGrid;

    template<Scalar T, bool isSpinPolarized>
    class BandGrid {
        using KPointGrid = ArrayND<KPoint<T, 0, isSpinPolarized>, 3>;
        using LatticeMatrix = CrystalCell<T>::LatticeMatrix;

        KPointGrid kPointGrid;
        LatticeMatrix repLatt;
        size_t electronCount;
        size_t numBand;
    public:
        BandGrid(LatticeMatrix repLatt_, size_t kPointX, size_t kPointY, size_t kPointZ, size_t electronCount_, size_t numBand_);
        BandGrid(const BandGrid&) = default;
        BandGrid(BandGrid&&) noexcept = default;
        ~BandGrid() = default;
        /* Operators */
        BandGrid& operator=(BandGrid band) noexcept { swap(band); return *this; }
        /* Operations */
        void swap(BandGrid& __restrict band) noexcept;
        /* Getters */
        [[nodiscard]] KPointGrid& getKPointGrid() { return kPointGrid; }
        [[nodiscard]] const KPointGrid& getKPointGrid() const noexcept { return kPointGrid; }
        [[nodiscard]] size_t getNumBand() const noexcept { return numBand; }
        [[nodiscard]] T getTotalEnergy() const noexcept;
        template<Vector V>
        [[nodiscard]] V getDensityOfStates(const V& atEnergy) const;
    private:
        [[nodiscard]] Vector3D<T> gradEnergy(size_t kPointId) const;
    };

    template<Scalar T, bool isSpinPolarized>
    BandGrid<T, isSpinPolarized>::BandGrid(LatticeMatrix repLatt_,
                                                    size_t kPointX,
                                                    size_t kPointY,
                                                    size_t kPointZ,
                                                    size_t electronCount_,
                                                    size_t numBand_)
            : kPointGrid({kPointX, kPointY, kPointZ})
            , repLatt(std::move(repLatt_))
            , electronCount(electronCount_)
            , numBand(numBand_) {
        assert(kPointGrid.getSize() != 0);
        assert(numBand >= (electronCount + 1) / 2);
        size_t kPointID = 0;

        const T kPointWeight = reciprocal(T(kPointGrid.getSize()));
        const T stepX = reciprocal(T(kPointX));
        const T stepY = reciprocal(T(kPointY));
        const T stepZ = reciprocal(T(kPointZ));

        Vector3D<T> k{};
        T& kx = k[0];
        T& ky = k[1];
        T& kz = k[2];

        kx = (T(1) - T(kPointX)) / T(2 * kPointX);
        for (size_t x = 1; x <= kPointX; ++x) {
            ky = (T(1) - T(kPointY)) / T(2 * kPointY);
            for (size_t y = 1; y <= kPointY; ++y) {
                kz = (T(1) - T(kPointZ)) / T(2 * kPointZ);
                for (size_t z = 1; z <= kPointZ; ++z) {
                    kPointGrid.asArray()[kPointID] = KPoint<T, 0, isSpinPolarized>(k, kPointWeight, numBand);
                    kz += stepZ;
                    ++kPointID;
                }
                ky += stepY;
            }
            kx += stepX;
        }
    }

    template<Scalar T, bool isSpinPolarized>
    void BandGrid<T, isSpinPolarized>::swap(BandGrid& __restrict band) noexcept {
        assert(this != &band && "[Error]: Self swap is likely a bug");
        swap(kPointGrid, band.kPointGrid);
        swap(electronCount, band.electronCount);
    }

    template<Scalar T, bool isSpinPolarized>
    T BandGrid<T, isSpinPolarized>::getTotalEnergy() const noexcept {
        T energy = T(0);
        for (const auto& kPoint : kPointGrid) {
            const T energyUp = kPoint.getBandEnergy(SpinState::Up).sum();
            const T energyDown = kPoint.getBandEnergy(SpinState::Down).head(electronCount / 2).sum();
            energy += (energyUp + energyDown) * kPoint.getWeight();
        }
        return energy;
    }
    /**
     * Reference:
     * [1] phys. stat. sol. (b) 179, 429-439 (2010); https://doi.org/10.1002/pssb.2221790218
     */
    template<Scalar T, bool isSpinPolarized>
    template<Vector V>
    V BandGrid<T, isSpinPolarized>::getDensityOfStates(const V& atEnergy) const {
        auto dos = V(atEnergy.getLength());
        for (size_t i = 0; i < atEnergy.getLength(); ++i) {
            const T energy = atEnergy[i];
            T density = T(0);
            for (size_t kPointId = 0; kPointId < kPointGrid.getSize(); ++i) {
                const T energy0 = kPointGrid[kPointId].getTotalEnergy();
                const auto gradE = gradEnergy(kPointId);
                const T normalizer = gradE[0] * gradE[1] * gradE[2] * T(0.5);

                T deltaDensity = T(0);
                for (int sigma1; sigma1 < 2; ++sigma1) {
                    for (int sigma2; sigma2 < 2; ++sigma2) {
                        for (int sigma3; sigma3 < 2; ++sigma3) {
                            const int sum = sigma1 + sigma2 + sigma3;
                            T temp = energy - energy0;
                            temp -= T(sigma1 == 0 ? 1 : -1) * gradE(0);
                            temp -= T(sigma2 == 0 ? 1 : -2) * gradE(1);
                            temp -= T(sigma3 == 0 ? 1 : -3) * gradE(2);
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

    template<Scalar T, bool isSpinPolarized>
    inline void swap(Physica::BandGrid<T, isSpinPolarized>& __restrict band1,
                     Physica::BandGrid<T, isSpinPolarized>& __restrict band2) noexcept {
        band1.swap(band2);
    }
    /**
     * \returns Gradient of energy in t-coordinate defined in [1]
     * 
     * Reference:
     * [1] phys. stat. sol. (b) 179, 429-439 (2010); https://doi.org/10.1002/pssb.2221790218
     */
    template<Scalar T, bool isSpinPolarized>
    Vector3D<T> BandGrid<T, isSpinPolarized>::gradEnergy(size_t kPointId) const {
        auto dimAdd = [](size_t dim, size_t dim_all) { return dim == dim_all - 1 ? 0 : dim + 1; };
        auto dimSub = [](size_t dim, size_t dim_all) { return dim == 0 ? dim_all - 1 : dim - 1; };
        const auto[x, y, z] = kPointGrid.indexToDim(kPointId);
        const size_t dimX = kPointGrid.getDimX();
        const size_t dimY = kPointGrid.getDimY();
        const size_t dimZ = kPointGrid.getDimZ();
        const T factor = T(0.25);
        const T gradX = (kPointGrid(dimAdd(x, dimX), y, z) -  kPointGrid(dimSub(x, dimX), y, z)) * factor;
        const T gradY = (kPointGrid(x, dimAdd(y, dimY), z) -  kPointGrid(x, dimSub(y, dimY), z)) * factor;
        const T gradZ = (kPointGrid(x, y, dimAdd(z, dimZ)) -  kPointGrid(x, y, dimSub(z, dimZ))) * factor;
        return {gradX, gradY, gradZ};
    }
}
