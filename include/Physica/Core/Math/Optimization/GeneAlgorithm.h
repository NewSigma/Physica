/*
 * Copyright 2020-2021 WeiBo He.
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

#include "Physica/Utils/Container/Array/Array.h"
#include "Physica/Core/MultiPrecision/Scalar.h"

namespace Physica::Core {
    /**
     * Search for minimum value
     */
    template<class Function, class VectorType>
    class GeneAlgorithm {
        using ScalarType = typename VectorType::ScalarType;
    public:
        struct AlgorithmConfig {
            float crossoverRate;
            float mutationRate;
            unsigned int populationSize;
            unsigned int maxGeneration;
            VectorType lowerBound;
            VectorType upperBound;
        };
    private:
        AlgorithmConfig config;
        Function func;
        Utils::Array<VectorType> population;
        Utils::Array<ScalarType> fitness;
        size_t minFitnessIndex;
    public:
        GeneAlgorithm(Function func_, const AlgorithmConfig& config_);
        ~GeneAlgorithm() = default;

        void solve();
        /* Getters */
        [[nodiscard]] const AlgorithmConfig& getConfig() const noexcept { return config; }
        [[nodiscard]] size_t getDim() const noexcept { return config.lowerBound.getLength(); }
        [[nodiscard]] const VectorType getOptimizedParams() const noexcept { return population[minFitnessIndex]; }
        [[nodiscard]] ScalarType getOptimizedValue() const noexcept { return fitness[minFitnessIndex]; }
        /* Setters */
        void setConfig(const AlgorithmConfig& c) { config = c; }
    private:
        void crossover();
        void mutation();
    };
}

#include "GeneAlgorithmImpl.h"
