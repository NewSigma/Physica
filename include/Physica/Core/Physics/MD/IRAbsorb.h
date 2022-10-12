/*
 * Copyright 2022 WeiBo He.
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

namespace Physica::Core {
    /**
     * Reference:
     * [1] S. Habershon, G. S. Fanourgakis, and D. E. Manolopoulos, J. Chem. Phys. 129, 074501 (2008).
     */
    template<class ScalarType>
    class IRAbsorb {
        using VectorType = Vector<ScalarType>;
    public:
        VectorType symmCorr;
        FFT<ScalarType, 1> fft;
        ScalarType factor;
        SavitzkyGolay<ScalarType> filter;
    public:
        IRAbsorb(const VectorType& dipoleCorr,
                 ScalarType deltaT,
                 double volume,
                 unsigned char filterRange,
                 size_t filterOrder);
        /* Operations */
        [[nodiscard]] VectorType makeWaveNum() const;
        [[nodiscard]] VectorType makeSpectrum() const;
        /* Getters */
        [[nodiscard]] size_t getDataLength() const noexcept { return symmCorr.getLength() / 2; }
        [[nodiscard]] ScalarType getDeltaWaveNum() const { return fft.getKSpaceDelta() / ScalarType(PhyConst<AU>::speedOfLight); }
    };

    template<class ScalarType>
    IRAbsorb<ScalarType>::IRAbsorb(const VectorType& dipoleCorr,
                                   ScalarType deltaT,
                                   double volume,
                                   unsigned char filterRange,
                                   size_t filterOrder)
            : symmCorr(2 * dipoleCorr.getLength())
            , fft(2 * dipoleCorr.getLength(), deltaT)
            , factor(M_PI / (3 * PhyConst<AU>::vacuumDielectric * PhyConst<AU>::speedOfLight * volume * PhyConst<AU>::boltzmannK * PhyConst<AU>::kToTemperature(298)))
            , filter(filterRange, filterRange, filterOrder, deltaT) {
        const Vector<ScalarType> temp = dipoleCorr.reverse();
        auto head = symmCorr.head(dipoleCorr.getLength());
        head = dipoleCorr;
        auto tail = symmCorr.tail(dipoleCorr.getLength());
        tail = temp;
        fft.transform(symmCorr);
    }

    template<class ScalarType>
    typename IRAbsorb<ScalarType>::VectorType IRAbsorb<ScalarType>::makeWaveNum() const {
        return VectorType::linspace(0, getDeltaWaveNum() * (getDataLength() - 1), getDataLength());
    }

    template<class ScalarType>
    typename IRAbsorb<ScalarType>::VectorType IRAbsorb<ScalarType>::makeSpectrum() const {
        const Vector<ScalarType> norm = toRealVector(fft.getKSpace()) * ScalarType(1 / (2 * M_PI));
        const ScalarType step =  getDeltaWaveNum();

        Vector<ScalarType> spectrum(getDataLength() + filter.getLRange() + filter.getRRange());
        size_t index = 0;
        for (size_t i = 0; i < filter.getLRange(); ++i, ++index)
            spectrum[index] = 0.0;
        for (size_t i = 0; i < getDataLength(); ++i, ++index)
            spectrum[index] = factor * abs(norm[i]) * square(step * i);
        for (size_t i = 0; i < filter.getRRange(); ++i, ++index)
            spectrum[index] = 0.0;
        filter.smooth(spectrum);
        return spectrum.segment(filter.getLRange(), getDataLength() + filter.getLRange());
    }
}
