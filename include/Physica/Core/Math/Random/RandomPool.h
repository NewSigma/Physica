/*
 * Copyright 2023-2024 WeiBo He.
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
    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed = Physica::Utils::Dynamic>
    class RandomPool;

    namespace Internal {
        template<class T> class Traits;

        template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
        class Traits<RandomPool<RandomGenerator, FixedSeed>> {
        public:
            constexpr static bool IsSeedFixed = FixedSeed != Utils::Dynamic;
        };
    }
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
    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
    class PHYSICA_API RandomPool {
        using This = RandomPool<RandomGenerator, FixedSeed>;
        using TraitsType = Internal::Traits<This>;
    public:
        using SeedType = typename RandomGenerator::result_type;
        using GeneratorType = RandomGenerator;
        constexpr static bool IsSeedFixed = TraitsType::IsSeedFixed;
    private:
        RandomGenerator gen;
        SeedType seed;
    public:
        /* Getters */
        [[nodiscard]] RandomGenerator& getGen() noexcept { return gen; }
        [[nodiscard]] inline SeedType getThreadSpecificSeed() const noexcept;
        /* Static members */
        [[nodiscard]] static RandomPool& getInstance();
    private:
        RandomPool();
        RandomPool(const RandomPool&) = default;
        RandomPool(RandomPool&&) noexcept = default;
        ~RandomPool() = default;
        /* Operators */
        RandomPool& operator=(const RandomPool&) = default;
        RandomPool& operator=(RandomPool&&) noexcept = default;
    };

    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
    RandomPool<RandomGenerator, FixedSeed>::RandomPool() {
        if constexpr (IsSeedFixed)
            seed = FixedSeed;
        else
            RandomSeed::rdrand(seed);
        gen.seed(getThreadSpecificSeed());
    }

    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
    inline typename RandomPool<RandomGenerator, FixedSeed>::SeedType
    RandomPool<RandomGenerator, FixedSeed>::getThreadSpecificSeed() const noexcept {
        if constexpr (IsSeedFixed) {
            const auto threadId = ThreadPool::getThreadInfo().id;
            return seed + threadId;
        }
        else
            return seed;
    }

    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
    RandomPool<RandomGenerator, FixedSeed>& RandomPool<RandomGenerator, FixedSeed>::getInstance() {
        thread_local static This instance{};
        return instance;
    }
}
