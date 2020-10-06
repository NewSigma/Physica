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
#ifndef PHYSICA_GENEALGORITHMIMPL_H
#define PHYSICA_GENEALGORITHMIMPL_H

namespace Physica::Core {
    template<ScalarType type>
    GeneAlgorithm<1, type>::GeneAlgorithm(
            const VectorFunction<type, false> &func,
            const Range &range,
            const AlgorithmConfig& config)
            : AbstractGeneAlgorithm(config)
            , func(func)
            , range(range)
            , regionLength(range.upperBound - range.lowerBound) {
        const auto population = config.population;
        const auto mode = config.mode;
        const auto& lower = range.lowerBound;
        if (mode == LinearChoose) {
            Scalar<type, false> stepLength =
                    regionLength / MultiScalar(static_cast<SignedScalarUnit>(population));
            Scalar<type, false> temp(lower);
            for (int i = 0; i < population; i++) {
                points[i] = temp;
                temp += stepLength;
            }
        }
        else {
            for (int i = 0; i < population; i++)
                points[i] = lower + randomScalar<MultiPrecision, false>() * regionLength;
        }
    }

    template<ScalarType type>
    Point<2, type> GeneAlgorithm<1, type>::solve() const {
        Q_UNUSED(type)
        unsigned int generation = 0;
        while (generation < config.maxGeneration) {
            crossover();
            mutation();
            generation += 1;
        }
        return points;
    }

    template<ScalarType type>
    void GeneAlgorithm<1, type>::crossover() {
        const auto population = config.population;
        const auto crossoverRate = config.crossoverRate;
        for (int i = 0; i < population; i++) {
            auto r = randomScalar<type, false>();
            if (double(r) < crossoverRate) {
                unsigned int randomIndex1 = random() % population;
                unsigned int randomIndex2 = random() % population;
                auto& random1 = points[randomIndex1];
                auto& random2 = points[randomIndex2];

                r = randomScalar<type, false>();
                //Whether random2 > random1 or not is not important.
                auto child = random1 + (random2 - random1) * r;
                auto child_y = func(child);
                auto y_random1 = func(random1);
                auto y_random2 = func(random2);
                if (child_y > y_random1)
                    random1 = child;
                else if (child_y > y_random2)
                    random2 = child;
            }
        }
    }

    template<ScalarType type>
    void GeneAlgorithm<1, type>::mutation() {
        Q_UNUSED(type)
        if(double(Physica::Core::randomScalar<MultiPrecision, false>()) < config.mutationRate) {
            unsigned int randomIndex = random() % config.population;
            points[randomIndex] =
                    range.lowerBound + randomScalar<MultiPrecision, false>() * regionLength;
        }
    }
}

#endif