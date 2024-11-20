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
#include "Physica/Core/Exception/MKL/VSL.h"
#include "Physica/Core/Parallel/ThreadPool.h"
#include "RandomSeed.h"

namespace Physica::Core {
    /**
     * \class Random provides a general, per-thread, reusable random generator implementation.
     * 
     * \tparam FixedSeed:
     * Multi-threads will break reproducibility and fixing random seed is not enough.
     * So \tparam FixedSeed is declared as a template param,
     * making it possible for other parts of the program to check whether the seed is fixed at compiling time.
     */
    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed = Physica::Dynamic>
    class PHYSICA_API Random {
        using This = Random<RandomGenerator, FixedSeed>;
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
        ~Random();
        /* Operators */
        [[nodiscard]] result_type operator()() { return gen(); }
        [[nodiscard]] operator VSLStreamStatePtr() noexcept { return pStream; }
        /* Getters */
        [[nodiscard]] RandomGenerator& getGen() noexcept { return gen; }
        [[nodiscard]] VSLStreamStatePtr getStream() noexcept { return pStream; }
        [[nodiscard]] inline SeedType getThreadSeed() const noexcept;
        /* Static members */
        [[nodiscard]] constexpr static result_type min() { return RandomGenerator::min(); }
        [[nodiscard]] constexpr static result_type max() { return RandomGenerator::max(); }

        [[nodiscard]] static Random& getInstance();
        [[nodiscard]] static Array<int> random_int(size_t length, int from, int to);
    private:
        Random();
        Random(const Random&) = default;
        Random(Random&&) noexcept = default;
        /* Operators */
        Random& operator=(const Random&) = default;
        Random& operator=(Random&&) noexcept = default;
        /* Static members */
        constexpr static int getMKLRngID();
    };

    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
    Random<RandomGenerator, FixedSeed>::Random() {
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
    Random<RandomGenerator, FixedSeed>::~Random() {
        vslCheck(vslDeleteStream(&pStream));
    }

    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
    inline typename Random<RandomGenerator, FixedSeed>::SeedType
    Random<RandomGenerator, FixedSeed>::getThreadSeed() const noexcept {
        if constexpr (IsSeedFixed) {
            const auto threadId = ThreadPool::getThreadInfo().id;
            return seed + threadId;
        }
        else
            return seed;
    }

    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
    Random<RandomGenerator, FixedSeed>& Random<RandomGenerator, FixedSeed>::getInstance() {
        thread_local static This instance{};
        return instance;
    }

    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
    Array<int> Random<RandomGenerator, FixedSeed>::random_int(size_t length, int from, int to) {
        assert(from <= to && to < INT_MAX);
        Array<int> result(length);
        vslCheck(viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, getInstance(), length, result.data(), from, to + 1));
        return result;
    }

    template<class RandomGenerator, typename RandomGenerator::result_type FixedSeed>
    constexpr int Random<RandomGenerator, FixedSeed>::getMKLRngID() {
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
    class Traits<Core::Random<RandomGenerator, FixedSeed>> {
    public:
        constexpr static bool IsSeedFixed = FixedSeed != Dynamic;
    };
}
