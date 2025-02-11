/*
 * Copyright 2024-2025 Weibo He.
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
#include "NumCharacter.h"

namespace Physica {
    /**
     * Compute correlation using FFT method as introduced in [1]
     * 
     * Reference:
     * [1] Computer Simulation of Liquids (2nd edn); https://doi.org/10.1093/oso/9780198803195.001.0001
     */
    template<Scalar T>
    class Correlation {
        using VectorType = VectorND<T>;

        FFT<T, 1> fft;
        VectorType corr;
        T mean;
        size_t numSample;
        size_t step;
    public:
        Correlation() = default;
        Correlation(size_t numStep);
        Correlation(const Correlation&) = default;
        Correlation(Correlation&&) noexcept = default;
        ~Correlation() = default;
        /* Operators */
        Correlation& operator=(Correlation obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void sample(T data);
        void resample() noexcept { step = 0; }
        [[nodiscard]] VectorType makeCorr(bool removeDrift) const;
        void swap(Correlation& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumStep() const noexcept { return corr.getLength(); }
        [[nodiscard]] T getMean() const noexcept { return mean; }
        [[nodiscard]] size_t getNumSample() const noexcept { return numSample; }
    };

    template<Scalar T>
    Correlation<T>::Correlation(size_t numStep)
            : fft(numStep * 2, PlanFlag::Estimate)
            , corr(numStep, 0)
            , numSample(0)
            , step(0) {
        assert(numStep > 1 && "[Error]: Invalid step number");
    }

    template<Scalar T>
    void Correlation<T>::sample(T data) {
        const size_t numStep = getNumStep();
        auto& rSpace = fft.getRSpace();
        auto head = rSpace.head(numStep);
        head[step] = data;
        step += 1;
        
        const bool isDataEnough = step == numStep;
        if (isDataEnough) {
            auto tail = rSpace.tail(numStep);
            tail = T(0);

            fft.transform();
            auto& kSpace = fft.getKSpace();
            toNextMean(mean, numSample, kSpace[0].real() / T(numStep));
            kSpace = kSpace.squaredNorms();
            fft.invTransform();

            for (size_t i = 0; i < numStep; ++i)
                head[i] /= T(numStep - i);
            toNextMean(corr, numSample, head);
            numSample += 1;
            step = 0;
        }
    }

    template<Scalar T>
    Correlation<T>::VectorType Correlation<T>::makeCorr(bool removeDrift) const {
        VectorType result = corr;
        if (removeDrift)
            result -= square(mean);
        return result;
    }

    template<Scalar T>
    void Correlation<T>::swap(Correlation& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        fft.swap(obj.fft);
        corr.swap(obj.corr);
        mean.swap(obj.mean);
        std::swap(numSample, obj.numSample);
        std::swap(step, obj.step);
    }
}
