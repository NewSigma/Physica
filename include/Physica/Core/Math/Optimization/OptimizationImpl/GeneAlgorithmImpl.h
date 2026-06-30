/*
 * Copyright 2020-2026 Weibo He.
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

namespace Physica {
    template<class Function, Vector V>
    GeneAlgorithm<Function, V>::GeneAlgorithm(Function func_, const AlgorithmConfig& config_)
            : config(config_)
            , func(func_)
            , population(config_.populationSize)
            , fitness(config_.populationSize)
            , minFitnessIndex(0) {
        assert(config.lowerBound.getLength() == config.upperBound.getLength());
        assert(config.populationSize > 0);
        for (size_t i = 0; i < population.getLength(); ++i) {
            population[i] = V::randomVector(config.lowerBound, config.upperBound);
            fitness[i] = func(std::cref(population[i]));
        }
    }

    template<class Function, Vector V>
    template<RNG R>
    void GeneAlgorithm<Function, V>::solve() {
        for (size_t iteration = 0; iteration < config.maxGeneration; ++iteration) {
            crossover<R>();
            mutation<R>();
        }
        minFitnessIndex = fitness.argmin();
    }

    template<class Function, Vector V>
    template<RNG R>
    void GeneAlgorithm<Function, V>::crossover() {
        const auto populationSize = config.populationSize;
        const auto crossoverRate = config.crossoverRate;
        for (size_t i = 0; i < populationSize; i++) {
            const ScalarType random = ScalarType::template random_uniform<R>();
            if (random < ScalarType(crossoverRate)) {
                unsigned int sampleIndex1 = std::rand() % populationSize;
                unsigned int sampleIndex2 = std::rand() % populationSize;
                auto& sample1 = population[sampleIndex1];
                auto& sample2 = population[sampleIndex2];

                const ScalarType random1 = ScalarType::template random_uniform<R>();
                //Whether random2 > random1 or not does not matter.
                V child = sample1 + (sample2 - sample1) * random1;
                ScalarType fitness_child = func(std::cref(child));
                if (fitness_child < fitness[sampleIndex1])
                    sample1 = std::move(child);
                else if (fitness_child < fitness[sampleIndex2])
                    sample2 = std::move(child);
            }
        }
    }

    template<class Function, Vector V>
    template<RNG R>
    void GeneAlgorithm<Function, V>::mutation() {
        const ScalarType random = ScalarType::template random_uniform<R>();
        if(random < ScalarType(config.mutationRate)) {
            unsigned int randomIndex = std::rand() % config.populationSize;
            population[randomIndex] = V::randomVector(config.lowerBound, config.upperBound);
        }
    }
}
