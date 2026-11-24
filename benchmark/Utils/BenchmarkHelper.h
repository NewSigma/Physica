/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Scalar/Real.h"
#include "Physica/Core/Math/Calculus/SpecialFunctions.h"
#include "Timer.h"

namespace Physica {
    class Benchmark {
        using ScalarType = float64;
    public:
        [[nodiscard]] static std::pair<ScalarType, ScalarType> run(std::invocable auto fn, unsigned int numTest, unsigned int numSample) {
            ScalarType mean = 0, var = 0;
            HostTimer timer;
            for (size_t i = 0; i < numTest; ++i) {
                ScalarType temp = 0;
                for (size_t j = 0; j < numSample; ++j) {
                    timer.tic();
                    fn();
                    temp.toNextMean(j, ScalarType(timer.toc()));
                }
                var.toNextVariance(mean, i, temp);
            }
            return std::make_pair(mean, sqrt(var));
        }
        /**
         * Suppose time1 ~ N(mean1, devia1^2) and time2 ~ N(mean2, devia2^2)
         */
        static ScalarType probabilityLarger1(ScalarType mean1, ScalarType devia1, ScalarType mean2, ScalarType devia2) {
            ScalarType mean = mean2 - mean1;
            ScalarType devia = sqrt(square(devia1) + square(devia2));
            return erfc(mean / (devia * sqrt(ScalarType(2.0)))) * 0.5;
        }
    };
}
