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
     * \class RandomPool Provides per-thread, reusable random generator support.
     * 
     * \tparam FixedSeed:
     * Multi-threads will break reproducibility and fixing random seed is not enough.
     * So \tparam FixedSeed is declared as a template param,
     * making it possible for other parts of the program to check whether the seed is fixed at compiling time.
     * 
     * Note:
     * If you use random numbers in a single thread code, construct a generator by hand is prefered.
     */
    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed = Physica::Utils::Dynamic>
    class PHYSICA_API RandomPool {
    public:
        using SeedType = typename RandomGenerator::result_type;
        using GeneratorType = RandomGenerator;
    public:
        /* Getters */
        [[nodiscard]] constexpr static bool isSeedFixed() { return FixedSeed != Utils::Dynamic; }
        [[nodiscard]] static RandomGenerator& getGen();
        [[nodiscard]] static SeedType getThreadSpecificSeed();
    };

    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
    RandomGenerator& RandomPool<RandomGenerator, FixedSeed>::getGen() {
        thread_local static RandomGenerator gen = RandomGenerator(getThreadSpecificSeed());
        return gen;
    }

    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
    typename RandomPool<RandomGenerator, FixedSeed>::SeedType RandomPool<RandomGenerator, FixedSeed>::getThreadSpecificSeed() {
        thread_local static SeedType seed = FixedSeed;
        if (isSeedFixed()) {
            const auto threadId = ThreadPool::getThreadInfo().id;
            return seed + threadId;
        }
        else {
            const bool isUninitialized = seed == Utils::Dynamic;
            if (isUninitialized)
                RandomSeed::rdrand(seed);
            return seed;
        }
    }
}
