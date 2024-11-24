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
    enum RandomOption {
        MT19937
    };

    class RandomBase {};

    template<class T>
    concept RandomGenerator = std::derived_from<T, RandomBase>;
    /**
     * \class Random provides a general, per-thread, reusable random generator implementation.
     * 
     * \tparam FixedSeed:
     * Multi-threads will break reproducibility and fixing random seed is not enough.
     * So \tparam FixedSeed is declared as a template param,
     * making it possible for other parts of the program to check whether the seed is fixed at compiling time.
     */
    template<RandomOption Option, uint64_t FixedSeed = Physica::Dynamic>
    class PHYSICA_API Random : public RandomBase {
        using This = Random<Option, FixedSeed>;
    public:
        using result_type = uint64_t;
        using SeedType = result_type;
        constexpr static bool IsSeedFixed = Traits<This>::IsSeedFixed;
    private:
        using GenType = std::mt19937;
    #ifndef PHYSICA_MKL
        using VSLStreamStatePtr = void*;
    #endif

        GenType gen;
        VSLStreamStatePtr pStream;

        SeedType seed;
    public:
        ~Random();
        /* Operators */
        [[nodiscard]] result_type operator()() { return gen(); }
        [[nodiscard]] operator VSLStreamStatePtr() noexcept { return pStream; }
        /* Operations */
        [[nodiscard("[Info]: Record the new seed for reproducible result")]] SeedType reseed();
        /* Getters */
        [[nodiscard]] GenType& getGen() noexcept { return gen; }
        [[nodiscard]] VSLStreamStatePtr getStream() noexcept { return pStream; }
        [[nodiscard]] SeedType getSeed() const noexcept { return seed; }
        /* Static members */
        [[nodiscard]] constexpr static result_type min() { return GenType::min(); }
        [[nodiscard]] constexpr static result_type max() { return GenType::max(); }

        [[nodiscard]] static Random& getInstance();
        [[nodiscard]] static Array<int> random_int(size_t length, int from, int to);
    private:
        Random();
        Random(const Random&) = default;
        Random(Random&&) noexcept = default;
        /* Operators */
        Random& operator=(const Random&) = default;
        Random& operator=(Random&&) noexcept = default;
        /* Getters */
        [[nodiscard]] inline SeedType getThreadSeed() const noexcept;
        /* Static members */
        constexpr static int getMKLRngID();
    };

    template<RandomOption Option, uint64_t FixedSeed>
    Random<Option, FixedSeed>::Random() {
        std::ignore = reseed();
    }

    template<RandomOption Option, uint64_t FixedSeed>
    Random<Option, FixedSeed>::~Random() {
        vslCheck(vslDeleteStream(&pStream));
    }

    template<RandomOption Option, uint64_t FixedSeed>
    typename Random<Option, FixedSeed>::SeedType Random<Option, FixedSeed>::reseed() {
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
        return seed;
    }

    template<RandomOption Option, uint64_t FixedSeed>
    Random<Option, FixedSeed>& Random<Option, FixedSeed>::getInstance() {
        thread_local static This instance{};
        return instance;
    }

    template<RandomOption Option, uint64_t FixedSeed>
    Array<int> Random<Option, FixedSeed>::random_int(size_t length, int from, int to) {
        assert(from <= to && to < INT_MAX);
        Array<int> result(length);
        vslCheck(viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, getInstance(), length, result.data(), from, to + 1));
        return result;
    }

    template<RandomOption Option, uint64_t FixedSeed>
    inline Random<Option, FixedSeed>::SeedType
    Random<Option, FixedSeed>::getThreadSeed() const noexcept {
        if constexpr (IsSeedFixed)
            return ThreadPool::isMainThread() ? seed : (seed + ThreadPool::getThreadInfo().id + 1);
        else
            return seed;
    }

    template<RandomOption Option, uint64_t FixedSeed>
    constexpr int Random<Option, FixedSeed>::getMKLRngID() {
    #ifdef PHYSICA_MKL
        return VSL_BRNG_MT19937;
    #else
        return 0;
    #endif
    }
}

namespace Physica {
    template<RandomOption Option, uint64_t FixedSeed>
    class Traits<Random<Option, FixedSeed>> {
    public:
        constexpr static bool IsSeedFixed = FixedSeed != Dynamic;
    };
}
