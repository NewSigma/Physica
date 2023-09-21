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
        mutable FFT<ScalarType, 1> fft;
        ScalarType factor;
        SavitzkyGolay<ScalarType> filter;
    public:
        IRAbsorb(VectorType dipoleCorr_,
                 ScalarType temperatureT,
                 ScalarType deltaT,
                 double volume,
                 unsigned char filterRange,
                 size_t filterOrder);
        IRAbsorb(const IRAbsorb&) = default;
        IRAbsorb(IRAbsorb&&) noexcept = default;
        ~IRAbsorb() = default;
        /* Operators */
        IRAbsorb& operator=(IRAbsorb obj) noexcept;
        /* Operations */
        [[nodiscard]] VectorType makeWaveNum() const;
        [[nodiscard]] VectorType makeSpectrum() const;
        void swap(IRAbsorb& obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getDataLength() const noexcept { return dipoleCorr.getLength(); }
        [[nodiscard]] ScalarType getDeltaWaveNum() const { return fft.getKSpaceDelta() / ScalarType(PhyConst<AU>::speedOfLight); }
    };

    template<class ScalarType>
    IRAbsorb<ScalarType>::IRAbsorb(VectorType dipoleCorr_,
                                   ScalarType temperatureT,
                                   ScalarType deltaT,
                                   double volume,
                                   unsigned char filterRange,
                                   size_t filterOrder)
            : dipoleCorr(std::move(dipoleCorr_))
            , fft(2 * dipoleCorr.getLength(), deltaT, PlanFlag::Estimate)
            , factor(M_PI / (3 * PhyConst<AU>::vacuumDielectric * PhyConst<AU>::speedOfLight * volume * PhyConst<AU>::boltzmannK * double(temperatureT)))
            , filter(filterRange, filterRange, filterOrder, deltaT) {}

    template<class ScalarType>
    IRAbsorb<ScalarType>& IRAbsorb<ScalarType>::operator=(IRAbsorb<ScalarType> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType>
    typename IRAbsorb<ScalarType>::VectorType IRAbsorb<ScalarType>::makeWaveNum() const {
        return VectorType::linspace(0, getDeltaWaveNum() * (getDataLength() - 1), getDataLength());
    }

    template<class ScalarType>
    typename IRAbsorb<ScalarType>::VectorType IRAbsorb<ScalarType>::makeSpectrum() const {
        Vector<ScalarType> symmCorr(2 * getDataLength());
        const Vector<ScalarType> temp = dipoleCorr.reverse();
        auto head = symmCorr.head(dipoleCorr.getLength());
        head = dipoleCorr;
        auto tail = symmCorr.tail(dipoleCorr.getLength());
        tail = temp;
        fft.transform(symmCorr);
        const Vector<ScalarType> norm = toRealVector(fft.getKSpace()) * ScalarType(1 / (2 * M_PI));
        const ScalarType step = fft.getKSpaceDelta();

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

    template<class ScalarType>
    void IRAbsorb<ScalarType>::swap(IRAbsorb& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        dipoleCorr.swap(obj.dipoleCorr);
        fft.swap(obj.fft);
        factor.swap(obj.factor);
        filter.swap(obj.filter);
    }
}
