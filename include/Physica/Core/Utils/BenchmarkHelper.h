/*
 * Copyright 2023-2024 Weibo He.
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

#include "Physica/Core/MultiPrecision/Real.h"
#include "Physica/Core/Math/Statistics/NumCharacter.h"
#include "Physica/Core/Math/Calculus/SpetialFunctions.h"
#include "Cycler.h"

namespace Physica::Core {
    class Benchmark {
        using ScalarType = float64;
    public:
        template<class Functor>
        [[nodiscard]] static std::pair<ScalarType, ScalarType> run(Functor func, unsigned int numTest, unsigned int numSample) {
            ScalarType mean = 0, var = 0;
            for (size_t i = 0; i < numTest; ++i) {
                ScalarType temp = 0;
                for (size_t j = 0; j < numSample; ++j) {
                    const auto from = Cycler::tic();
                    func();
                    const auto to = Cycler::toc();
                    toNextMean(temp, j, ScalarType(Cycler::toSeconds(to - from)));
                }
                toNextVariance(var, mean, i, temp);
            }
            return std::make_pair(mean, sqrt(var));
        }
        /**
         * Suppose time1 ~ N(mean1, devia1^2) and time2 ~ N(mean2, devia2^2)
         */
        static ScalarType probabilityLarger1(ScalarType mean1, ScalarType devia1, ScalarType mean2, ScalarType devia2) {
            ScalarType mean = mean2 - mean1;
            ScalarType devia = sqrt(square(devia1) + square(devia2));
            return Core::erfc(mean / (devia * sqrt(ScalarType(2.0)))) * 0.5;
        }
    };
}
