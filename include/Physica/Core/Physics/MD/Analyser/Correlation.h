/*
 * Copyright 2024 WeiBo He.
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

#include "Physica/Core/Math/Transform/FFT.h"

namespace Physica::Core {
    /**
     * Compute correlation using FFT method as introduced in [1]
     * 
     * Reference:
     * [1] Computer Simulation of Liquids (2nd edn); https://doi.org/10.1093/oso/9780198803195.001.0001
     */
    template<class ScalarType>
    class Correlation {
        using VectorType = Vector<ScalarType>;

        FFT<ScalarType, 1> fft;
        VectorType corr;
        size_t numSample;
    public:
        Correlation(size_t numStep);
        Correlation(const Correlation&) = default;
        Correlation(Correlation&&) noexcept = default;
        ~Correlation() = default;
        /* Operators */
        Correlation& operator=(Correlation obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void sample(const VectorType& v);
        void swap(Correlation& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumStep() const noexcept { return corr.getLength(); }
        [[nodiscard]] const VectorType& getCorr() const noexcept { return corr; }
        [[nodiscard]] size_t getNumSample() const noexcept { return numSample; }
    };

    template<class ScalarType>
    Correlation<ScalarType>::Correlation(size_t numStep)
            : fft(numStep * 2, 1, PlanFlag::Estimate)
            , corr(numStep, 0)
            , numSample(0) {}

    template<class ScalarType>
    void Correlation<ScalarType>::sample(const VectorType& v) {
        const size_t numStep = getNumStep();
        auto& rSpace = fft.getRSpace();
        auto head = rSpace.head(numStep);
        auto tail = rSpace.tail(numStep);
        head = v;
        tail = ScalarType(0);

        fft.transform();
        auto& kSpace = fft.getKSpace();
        kSpace = toNormVector(kSpace);
        fft.invTransform();

        const ScalarType factor = reciprocal(ScalarType(fft.getRSpaceSize()));
        auto head = rSpace.head(numStep);
        head *= factor;
        toNextMean(corr, numSample, head);
        numSample += 1;
    }

    template<class ScalarType>
    void Correlation<ScalarType>::swap(Correlation& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        fft.swap(obj.fft);
        corr.swap(obj.corr);
    }
}
