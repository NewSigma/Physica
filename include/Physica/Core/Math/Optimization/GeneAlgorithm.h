/*
 * Copyright 2020 WeiBo He.
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
#ifndef PHYSICA_GENEALGORITHM_H
#define PHYSICA_GENEALGORITHM_H

#include <ctime>
#include "Physica/Core/Math/Calculus/Function/VectorFunction/VectorFunction.h"
#include "Physica/Core/Math/Geometry/Point.h"

namespace Physica::Core {
    class AbstractGeneAlgorithm {
    public:
        enum ChooseMode {
            LinearChoose,
            RandomChoose
        };

        struct AlgorithmConfig {
            float crossoverRate = 0.6;
            float mutationRate = 0.1;
            //The size of initial points we will spawn.
            unsigned int population = 100;
            //Stopping args
            unsigned int maxGeneration = 100;
            ChooseMode mode = LinearChoose;
        };
    protected:
        AlgorithmConfig config;
    public:
        ~AbstractGeneAlgorithm() = default;
        /* Getters */
        [[nodiscard]] const AlgorithmConfig& getConfig() const { return config; }
        /* Setters */
        void setConfig(const AlgorithmConfig& c) { config = c; }
    protected:
        explicit AbstractGeneAlgorithm(const AlgorithmConfig& c);
    };
    
    template<size_t dim, ScalarOption option>
    class GeneAlgorithm : public AbstractGeneAlgorithm {
        static_assert(dim > 0, "dim must be positive.");
    };
    
    template<ScalarOption option>
    class GeneAlgorithm<1, option> : public AbstractGeneAlgorithm {
    public:
        struct Range {
            Scalar<option, false> lowerBound;
            Scalar<option, false> upperBound;
        };
    private:
        VectorFunction<option, false> func;
        Range range;
        Scalar<option, false> regionLength;
        std::vector<Scalar<option, false>> points;
    public:
        GeneAlgorithm(const VectorFunction<option, false>& func, const Range& range, const AlgorithmConfig& config);
        ~GeneAlgorithm() = default;

        Point<2, option> solve() const;
    private:
        void crossover();
        void mutation();
    };
}

#include "GeneAlgorithmImpl.h"

#endif