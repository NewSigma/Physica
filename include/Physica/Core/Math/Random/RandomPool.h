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

#include "RandomSeed.h"
#include "Physica/Core/Parallel/ThreadPool.h"

namespace Physica::Core {
    /**
     * Provides per-thread, reusable random generator support.
     * 
     * If you do not use random numbers in a multi-thread code, construct a generator by hand is prefered.
     */
    template<class RandomGenerator>
    class RandomPool {
    public:
        using SeedType = typename RandomGenerator::result_type;
        using GeneratorType = RandomGenerator;

        static SeedType FixedSeed;
    public:
        /* Getters */
        [[nodiscard]] static RandomGenerator& getGen();
        [[nodiscard]] static bool isSeedFixed() { return FixedSeed != Utils::Dynamic; }
        [[nodiscard]] static SeedType getThreadSpecificSeed();
    private:
        static RandomGenerator makeGenerator();
    };

    template<class RandomGenerator>
    typename RandomPool<RandomGenerator>::SeedType RandomPool<RandomGenerator>::FixedSeed = Utils::Dynamic;

    template<class RandomGenerator>
    RandomGenerator& RandomPool<RandomGenerator>::getGen() {
        thread_local static RandomGenerator gen = makeGenerator();
        return gen;
    }

    template<class RandomGenerator>
    typename RandomPool<RandomGenerator>::SeedType RandomPool<RandomGenerator>::getThreadSpecificSeed() {
        if (!isSeedFixed())
            return Utils::Dynamic;
        const auto threadId = ThreadPool::getThreadInfo().id;
        return FixedSeed + threadId;
    }

    template<class RandomGenerator>
    RandomGenerator RandomPool<RandomGenerator>::makeGenerator() {
        SeedType seed;
        if (isSeedFixed())
            seed = getThreadSpecificSeed();
        else
            RandomSeed::rdrand(seed);
        return RandomGenerator(seed);
    }
}
