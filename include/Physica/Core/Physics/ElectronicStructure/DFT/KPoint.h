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

namespace Physica::Core {
    template<class ScalarType, size_t NumBand, bool isSpinPolarized>
    class KPoint {
    public:
        constexpr static size_t NumSpin = isSpinPolarized ? 2 : 1;
        using ComplexType = Complex<ScalarType>;
        using Vector3D = Vector3D<ScalarType>;
        using BandEnergy = Vector<ScalarType, NumBand>;
        using BandEnergyPair = std::pair<BandEnergy, BandEnergy>;
    private:
        Vector3D pos;
        ScalarType weight;
        Array<BandEnergy, NumSpin> bandE;
    public:
        KPoint() = default;
        KPoint(Vector3D pos_, ScalarType weight_, size_t numBand);
        KPoint(const KPoint&) = default;
        KPoint(KPoint&&) noexcept = default;
        ~KPoint() = default;
        /* Operators */
        KPoint& operator=(KPoint k) noexcept;
        /* Operations */
        void swap(KPoint& __restrict kPoint) noexcept;
        /* Getters */
        [[nodiscard]] const Vector3D& getPos() const noexcept { return pos; }
        [[nodiscard]] const ScalarType& getWeight() const noexcept { return weight; }
        [[nodiscard]] inline const BandEnergy& getBandEnergy(SpinState spin) const noexcept;
        /* Setters */
        template<class VectorType>
        void setBandEnergy(SpinState spin, const RValueVector<VectorType>& v);
    };

    template<class ScalarType, size_t NumBand, bool isSpinPolarized>
    KPoint<ScalarType, NumBand, isSpinPolarized>::KPoint(Vector3D pos_, ScalarType weight_, size_t numBand)
            : pos(std::move(pos_)), weight(std::move(weight_)), bandE(NumSpin, numBand) {}

    template<class ScalarType, size_t NumBand, bool isSpinPolarized>
    KPoint<ScalarType, NumBand, isSpinPolarized>& KPoint<ScalarType, NumBand, isSpinPolarized>::operator=(KPoint k) noexcept {
        swap(k);
        return *this;
    }

    template<class ScalarType, size_t NumBand, bool isSpinPolarized>
    void KPoint<ScalarType, NumBand, isSpinPolarized>::swap(KPoint& __restrict kPoint) noexcept {
        assert(this != &kPoint && "[Error]: Self swap is likely a bug");
        pos.swap(kPoint.pos);
        weight.swap(kPoint.weight);
        bandE.swap(kPoint.bandE);
    }

    template<class ScalarType, size_t NumBand, bool isSpinPolarized>
    inline const typename KPoint<ScalarType, NumBand, isSpinPolarized>::BandEnergy&
    KPoint<ScalarType, NumBand, isSpinPolarized>::getBandEnergy(SpinState spin) const noexcept {
        return bandE[isSpinPolarized ? int(spin) : 0];
    }

    template<class ScalarType, size_t NumBand, bool isSpinPolarized>
    template<class VectorType>
    void KPoint<ScalarType, NumBand, isSpinPolarized>::setBandEnergy(SpinState spin, const RValueVector<VectorType>& v) {
        BandEnergy& energy = bandE[isSpinPolarized ? int(spin) : 0];
        const size_t length = energy.getLength();
        for (size_t i = 0; i < length; ++i)
            energy[i] = v.calc(i);
    }
}
