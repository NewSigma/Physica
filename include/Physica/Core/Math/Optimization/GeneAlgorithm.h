/*
 * Copyright 2020-2024 Weibo He.
 *
 * This file is part of Physica.

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

#include "Physica/Core/Utils/Container/Array.h"
#include "Physica/Core/Scalar/Real.h"

namespace Physica::Core {
    /**
     * Search for minimum value
     */
    template<class Function, Vector T>
    class GeneAlgorithm {
        using ScalarType = T::ScalarType;
    public:
        struct AlgorithmConfig {
            float crossoverRate;
            float mutationRate;
            unsigned int populationSize;
            unsigned int maxGeneration;
            T lowerBound;
            T upperBound;
        };
    private:
        AlgorithmConfig config;
        Function func;
        Array<T> population;
        Array<ScalarType> fitness;
        size_t minFitnessIndex;
    public:
        GeneAlgorithm(Function func_, const AlgorithmConfig& config_);
        ~GeneAlgorithm() = default;

        template<RandomGenerator R>
        void solve();
        /* Getters */
        [[nodiscard]] const AlgorithmConfig& getConfig() const noexcept { return config; }
        [[nodiscard]] size_t getDim() const noexcept { return config.lowerBound.getLength(); }
        [[nodiscard]] const T getOptimizedParams() const noexcept { return population[minFitnessIndex]; }
        [[nodiscard]] ScalarType getOptimizedValue() const noexcept { return fitness[minFitnessIndex]; }
        /* Setters */
        void setConfig(const AlgorithmConfig& c) { config = c; }
    private:
        template<RandomGenerator R>
        void crossover();
        template<RandomGenerator R>
        void mutation();
    };
}

#include "OptimizationImpl/GeneAlgorithmImpl.h"
