/*
 * Copyright 2026 Weibo He.
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

#include "GreenSampler/ScalarSampler.h"

namespace Physica {
    /**
     * A direct implementation for adaptive doping
     */
    template<Scalar T, RNG R>
    void fastDoping(int iteration, int windowSize, T targetRho, T stepsize, HubbardParams<T>& params, auto& dqmc, auto&&... args) {
        assert(windowSize % 2 == 1 && "[Error]: Required to avoid divide by zero");
        assert(T(0) <= targetRho && targetRho <= T(2) && "[Error]: Invalid target");
        ScalarSampler<T> sampler1(params, windowSize);
        ScalarSampler<T> sampler2(params, windowSize);
        for (int i = 0; i < windowSize + iteration; ++i) {
            bool initialized = i >= windowSize;
            if (initialized) {
                T rho = sampler1.calcMean();
                assert(rho.isPositive());
                T compressibility = T(params.getBeta() * params.getNumSite()) * (sampler2.calcMean() - square(rho));
                T mu = params.getChemMu() + (T(targetRho) - rho) / compressibility * stepsize;
                params.setChemMu(mu);
            }

            dqmc.template step<R>(std::forward<decltype(args)>(args)...);
            sampler1.sample(dqmc.getGreens(), dqmc.getRSign(), ScalarSampler<T>::Density);
            sampler2.sample(dqmc.getGreens(), dqmc.getRSign(), ScalarSampler<T>::Density2);
        }
    }
}
