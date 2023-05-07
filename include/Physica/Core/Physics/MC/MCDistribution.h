/*
 * Copyright 2023 WeiBo He.
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

namespace Physica::Core {
    template<class ScalarType>
    class MCDistribution {
        using VectorType = Vector<ScalarType>;

        VectorType samples;
        std::uniform_int_distribution<size_t> dist;
    public:
        template<class Functor, class RandomGenerator>
        MCDistribution makeDistribution(Functor func, RandomGenerator& gen, size_t numSample);
    private:
        MCDistribution(VectorType samples_);
    };

    template<class ScalarType>
    MCDistribution<ScalarType>::MCDistribution(VectorType samples_) : samples(std::move(samples_)) {
        dist = std::uniform_int_distribution<size_t>(0, samples.getLength() - 1);
    }

    template<class ScalarType>
    template<class Functor, class RandomGenerator>
    MCDistribution<ScalarType> MCDistribution<ScalarType>::makeDistribution(Functor func, RandomGenerator& gen, size_t numSample) {
        VectorType samples(numSample);
        ScalarType x = 0;
        ScalarType p = 1;
        std::uniform_real_distribution<> step_dist(-0.5, 0.5);
        std::uniform_real_distribution<> uniform_dist{};
        for (size_t i = 0; i < samples.getLength(); ++i) {
            ScalarType new_x = x + step_dist(gen);
            ScalarType new_p = func(new_x);
            if (new_p > p) {
                x = new_x;
                p = new_p;
            }
            else {
                ScalarType temp = uniform_dist(gen);
                if (temp < new_p / p) {
                    x = new_x;
                    p = new_p;
                }
            }
            samples[i] = x;
        }
        return MCDistribution(std::move(samples));
    }
}
