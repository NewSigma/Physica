/*
 * Copyright 2022-2023 WeiBo He.
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

#include "Physica/Core/Physics/PhyConst.h"
#include "Physica/Core/Physics/Filter/SavitzkyGolay.h"
#include "Physica/Core/Math/Transform/FFT.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] S. Habershon, G. S. Fanourgakis, and D. E. Manolopoulos, J. Chem. Phys. 129, 074501 (2008).
     */
    template<class ScalarType>
    class IRAbsorb {
        static_assert(is_scalar<ScalarType>::value, "[Error]: ScalarType must be a scalar");
        using VectorType = Vector<ScalarType>;
    public:
        VectorType dipoleCorr;
        FFT<ScalarType, 1> fft;
        ScalarType factor;
        SavitzkyGolay<ScalarType> filter;
    public:
        IRAbsorb(VectorType dipoleCorr_,
                 ScalarType temperatureT,
                 ScalarType deltaT,
                 ScalarType volume,
                 unsigned char filterRange,
                 size_t filterOrder);
        IRAbsorb(const IRAbsorb&) = default;
        IRAbsorb(IRAbsorb&&) noexcept = default;
        ~IRAbsorb() = default;
        /* Operators */
        IRAbsorb& operator=(IRAbsorb obj) noexcept;
        /* Operations */
        [[nodiscard]] VectorType makeWaveNum() const;
        [[nodiscard]] VectorType makeSpectrum();
        void swap(IRAbsorb& obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getDataLength() const noexcept { return dipoleCorr.getLength(); }
        [[nodiscard]] ScalarType getDeltaOmegaW() const noexcept { return fft.getKSpaceDelta(); }
        [[nodiscard]] ScalarType getDeltaWaveNum() const noexcept { return getDeltaOmegaW() / ScalarType(2 * M_PI * PhyConst<AU>::speedOfLight); }
    };

    template<class ScalarType>
    IRAbsorb<ScalarType>::IRAbsorb(VectorType dipoleCorr_,
                                   ScalarType temperatureT,
                                   ScalarType deltaT,
                                   ScalarType volume,
                                   unsigned char filterRange,
                                   size_t filterOrder)
            : dipoleCorr(std::move(dipoleCorr_))
            , fft(2 * dipoleCorr.getLength() - 1, deltaT, PlanFlag::Estimate)
            , factor(1 / (6 * PhyConst<AU>::vacuumDielectric * PhyConst<AU>::speedOfLight * double(volume) * PhyConst<AU>::boltzmannK * double(temperatureT)))
            , filter(filterRange, filterRange, filterOrder, deltaT) {}

    template<class ScalarType>
    IRAbsorb<ScalarType>& IRAbsorb<ScalarType>::operator=(IRAbsorb<ScalarType> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType>
    typename IRAbsorb<ScalarType>::VectorType IRAbsorb<ScalarType>::makeWaveNum() const {
        const size_t size = fft.getKSpaceSize();
        return VectorType::linspace(0, getDeltaWaveNum() * (size - 1), size);
    }

    template<class ScalarType>
    typename IRAbsorb<ScalarType>::VectorType IRAbsorb<ScalarType>::makeSpectrum() {
        auto& rSpace = fft.getRSpace();
        rSpace.head(dipoleCorr.getLength()) = dipoleCorr;
        rSpace.tail(dipoleCorr.getLength()) = dipoleCorr.tail(1).reverse();
        fft.transform();

        const auto& kSpace = fft.getKSpace();
        const size_t kSpaceSize = fft.getKSpaceSize();
        const ScalarType deltaOmegaW = getDeltaOmegaW();
        Vector<ScalarType> spectrum(kSpaceSize + filter.getLRange() + filter.getRRange());
        size_t index = 0;
        for (size_t i = 0; i < filter.getLRange(); ++i, ++index)
            spectrum[index] = 0.0;
        for (size_t i = 0; i < kSpaceSize; ++i, ++index)
            spectrum[index] = factor * abs(kSpace[i].getReal()) * square(deltaOmegaW * i);
        for (size_t i = 0; i < filter.getRRange(); ++i, ++index)
            spectrum[index] = 0.0;
        filter.smooth(spectrum);
        return spectrum.segment(filter.getLRange(), kSpaceSize + filter.getLRange());
    }

    template<class ScalarType>
    void IRAbsorb<ScalarType>::swap(IRAbsorb& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        dipoleCorr.swap(obj.dipoleCorr);
        fft.swap(obj.fft);
        factor.swap(obj.factor);
        filter.swap(obj.filter);
    }
}
