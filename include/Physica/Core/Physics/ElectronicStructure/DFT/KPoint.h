/*
 * Copyright 2022-2024 Weibo He.
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

#include "SpinPair.h"

namespace Physica {
    template<Scalar T, size_t NumBand, bool isSpinPolarized>
    class KPoint {
    public:
        constexpr static size_t NumSpin = isSpinPolarized ? 2 : 1;
        using ComplexType = Complex<T>;
        using BandEnergy = DenseVector<T, NumBand>;
        using BandEnergyPair = std::pair<BandEnergy, BandEnergy>;
    private:
        Vector3D<T> pos;
        T weight;
        Array<BandEnergy, NumSpin> bandE;
    public:
        KPoint() = default;
        KPoint(Vector3D<T> pos_, T weight_, size_t numBand);
        KPoint(const KPoint&) = default;
        KPoint(KPoint&&) noexcept = default;
        ~KPoint() = default;
        /* Operators */
        KPoint& operator=(KPoint k) noexcept;
        /* Operations */
        void swap(KPoint& __restrict kPoint) noexcept;
        /* Getters */
        [[nodiscard]] const Vector3D<T>& getPos() const noexcept { return pos; }
        [[nodiscard]] const T& getWeight() const noexcept { return weight; }
        [[nodiscard]] const BandEnergy& getBandEnergy(SpinState spin) const noexcept;
        /* Setters */
        void setBandEnergy(SpinState spin, const Vector auto& v);
    };

    template<Scalar T, size_t NumBand, bool isSpinPolarized>
    KPoint<T, NumBand, isSpinPolarized>::KPoint(Vector3D<T> pos_, T weight_, size_t numBand)
            : pos(std::move(pos_)), weight(std::move(weight_)), bandE(NumSpin, numBand) {}

    template<Scalar T, size_t NumBand, bool isSpinPolarized>
    KPoint<T, NumBand, isSpinPolarized>& KPoint<T, NumBand, isSpinPolarized>::operator=(KPoint k) noexcept {
        swap(k);
        return *this;
    }

    template<Scalar T, size_t NumBand, bool isSpinPolarized>
    void KPoint<T, NumBand, isSpinPolarized>::swap(KPoint& __restrict kPoint) noexcept {
        assert(this != &kPoint && "[Error]: Self swap is likely a bug");
        pos.swap(kPoint.pos);
        weight.swap(kPoint.weight);
        bandE.swap(kPoint.bandE);
    }

    template<Scalar T, size_t NumBand, bool isSpinPolarized>
    const KPoint<T, NumBand, isSpinPolarized>::BandEnergy&
    KPoint<T, NumBand, isSpinPolarized>::getBandEnergy(SpinState spin) const noexcept {
        return bandE[isSpinPolarized ? int(spin) : 0];
    }

    template<Scalar T, size_t NumBand, bool isSpinPolarized>
    void KPoint<T, NumBand, isSpinPolarized>::setBandEnergy(SpinState spin, const Vector auto& v) {
        BandEnergy& energy = bandE[isSpinPolarized ? int(spin) : 0];
        const size_t length = energy.getLength();
        for (size_t i = 0; i < length; ++i)
            energy[i] = v.calc(i);
    }
}
