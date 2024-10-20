/*
 * Copyright 2023-2024 Weibo He.
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

#include <random>
#include <Physica/Core/Exception/MKL/VSL.h>
#include <Physica/Core/Parallel/ThreadPool.h>
#include "RandomSeed.h"

namespace Physica::Core {
    /**
     * \class RandomPool provides per-thread, reusable random generator support.
     * 
     * \tparam FixedSeed:
     * Multi-threads will break reproducibility and fixing random seed is not enough.
     * So \tparam FixedSeed is declared as a template param,
     * making it possible for other parts of the program to check whether the seed is fixed at compiling time.
     */
    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed = Physica::Dynamic>
    class PHYSICA_API RandomPool {
        using This = RandomPool<RandomGenerator, FixedSeed>;
    public:
        using result_type = typename RandomGenerator::result_type;
        using SeedType = result_type;
        using GeneratorType = RandomGenerator;
        constexpr static bool IsSeedFixed = Traits<This>::IsSeedFixed;
    private:
    #ifndef PHYSICA_MKL
        using VSLStreamStatePtr = void*;
    #endif

        RandomGenerator gen;
        VSLStreamStatePtr pStream;

        SeedType seed;
    public:
        ~RandomPool();
        /* Operators */
        [[nodiscard]] result_type operator()() { return gen(); }
        [[nodiscard]] operator VSLStreamStatePtr() noexcept { return pStream; }
        /* Getters */
        [[nodiscard]] RandomGenerator& getGen() noexcept { return gen; }
        [[nodiscard]] VSLStreamStatePtr getStream() noexcept { return pStream; }
        [[nodiscard]] inline SeedType getThreadSeed() const noexcept;
        /* Static members */
        [[nodiscard]] static RandomPool& getInstance();
    private:
        RandomPool();
        RandomPool(const RandomPool&) = default;
        RandomPool(RandomPool&&) noexcept = default;
        /* Operators */
        RandomPool& operator=(const RandomPool&) = default;
        RandomPool& operator=(RandomPool&&) noexcept = default;
        /* Static members */
        constexpr static int getMKLRngID();
    };

    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
    RandomPool<RandomGenerator, FixedSeed>::RandomPool() {
        if constexpr (IsSeedFixed)
            seed = FixedSeed;
        else
            RandomSeed::rdrand(seed);

        auto tseed = getThreadSeed();
        gen.seed(tseed);
        if constexpr (HasMKL()) {
            RandomSeed::toNextSeed(tseed);
            vslCheck(vslNewStream(&pStream, getMKLRngID(), tseed));
        }
    }

    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
    RandomPool<RandomGenerator, FixedSeed>::~RandomPool() {
        vslCheck(vslDeleteStream(&pStream));
    }

    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
    inline typename RandomPool<RandomGenerator, FixedSeed>::SeedType
    RandomPool<RandomGenerator, FixedSeed>::getThreadSeed() const noexcept {
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

    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
    constexpr int RandomPool<RandomGenerator, FixedSeed>::getMKLRngID() {
        static_assert(std::is_same<RandomGenerator, std::mt19937>::value, "[Error]: Other RandomGenerator not implemented");
    #ifdef PHYSICA_MKL
        return VSL_BRNG_MT19937;
    #else
        return 0;
    #endif
    }
}

namespace Physica {
    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
    class Traits<Core::RandomPool<RandomGenerator, FixedSeed>> {
    public:
        constexpr static bool IsSeedFixed = FixedSeed != Dynamic;
    };
}
