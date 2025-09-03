/*
 * Copyright 2023-2025 Weibo He.
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
#ifdef PHYSICA_CUDA
    #include "Physica/Core/Exception/CUDA/cuRAND.cuh"
#endif
#include "Physica/Core/Exception/MKL/VSL.h"
#include "Physica/Core/Parallel/ThreadPool.h"
#include "RandomSeed.h"

namespace Physica {
    namespace Internal {
        class RandomBase {};

        class QRandomBase : public RandomBase {};
    }

    enum RandomOption {
        MT19937
    };
    /**
     * Random Number Generator
     */
    template<class T>
    concept RNG = std::derived_from<T, Internal::RandomBase>;
    /**
     * Quasi Random Number Generator
     */
    template<class T>
    concept QRNG = std::derived_from<T, Internal::QRandomBase>;
    /**
     * \class Random provides a general, per-thread, reusable random generator implementation.
     * 
     * \tparam FixedSeed:
     * Multi-threads will break reproducibility and fixing random seed is not enough.
     * So \tparam FixedSeed is declared as a template param,
     * making it possible for other parts of the program to check whether the seed is fixed at compiling time.
     */
    template<RandomOption Option, uint64_t FixedSeed = Physica::Dynamic>
    class PHYSICA_API Random : public Internal::RandomBase {
        using This = Random<Option, FixedSeed>;
    public:
        using result_type = uint64_t;
        using SeedType = result_type;
        constexpr static bool IsSeedFixed = Traits<This>::IsSeedFixed;
    private:
        using GenType = std::mt19937;

        GenType gen;
        VSLStreamStatePtr pStream;
        [[no_unique_address]] curandGenerator_t curand;

        SeedType seed;
    public:
        ~Random();
        /* Operators */
        [[nodiscard]] result_type operator()() { return gen(); }
        [[nodiscard]] operator VSLStreamStatePtr() noexcept { return pStream; }
        [[nodiscard]] operator curandGenerator_t() noexcept { return curand; }
        /* Operations */
        [[nodiscard("[Info]: Record the new seed for reproducible result")]] SeedType reseed();
        void reseed(SeedType seed_);
        /* Getters */
        [[nodiscard]] GenType& getGen() noexcept { return gen; }
        [[nodiscard]] VSLStreamStatePtr getStream() noexcept { return pStream; }
        [[nodiscard]] SeedType getSeed() const noexcept { return seed; }
        /* Static members */
        [[nodiscard]] constexpr static result_type min() { return GenType::min(); }
        [[nodiscard]] constexpr static result_type max() { return GenType::max(); }

        [[nodiscard]] static This& getInstance() noexcept;
        [[nodiscard]] static Array<int> random_int(size_t length, int from, int to);
        static void random_int(Array<int>& arr, int from, int to);
    private:
        Random() noexcept;
        Random(const This&) = default;
        Random(This&&) noexcept = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Getters */
        [[nodiscard]] SeedType getThreadSeed() const noexcept;
        /* Static members */
        constexpr static int getMKLRngID();
        constexpr static curandRngType_t curandRngID();
    };

    template<RandomOption Option, uint64_t FixedSeed>
    Random<Option, FixedSeed>::Random() noexcept {
        try {
            std::ignore = reseed();
        }
        catch (...) {
            fprintf(stderr, "RNG init failed");
            std::abort();
        }
    }

    template<RandomOption Option, uint64_t FixedSeed>
    Random<Option, FixedSeed>::~Random() {
        check_vsl(vslDeleteStream(&pStream));
    #ifdef PHYSICA_CUDA
        check(curandDestroyGenerator(curand));
    #endif
    }

    template<RandomOption Option, uint64_t FixedSeed>
    auto Random<Option, FixedSeed>::reseed() -> SeedType {
        if constexpr (IsSeedFixed)
            reseed(FixedSeed);
        else {
            RandomSeed::rdrand(seed);
            reseed(seed);
        }
        return seed;
    }

    template<RandomOption Option, uint64_t FixedSeed>
    void Random<Option, FixedSeed>::reseed(SeedType seed_) {
        if constexpr (IsSeedFixed) {
            assert(seed_ == FixedSeed);
            seed = FixedSeed;
        }
        else
            seed = seed_;

        auto tseed = getThreadSeed();
        gen.seed(tseed);
        if constexpr (HasMKL()) {
            RandomSeed::toNextSeed(tseed);
            check_vsl(vslNewStream(&pStream, getMKLRngID(), tseed));
        }
        if constexpr (HasCUDA()) {
        #ifdef PHYSICA_CUDA
            RandomSeed::toNextSeed(tseed);
            check(curandCreateGenerator(&curand, curandRngID()));
            check(curandSetGeneratorOrdering(curand, CURAND_ORDERING_PSEUDO_DYNAMIC));
            check(curandSetPseudoRandomGeneratorSeed(curand, tseed));
        #endif
        }
    }

    template<RandomOption Option, uint64_t FixedSeed>
    auto Random<Option, FixedSeed>::getInstance() noexcept -> This& {
        thread_local static This instance{};
        return instance;
    }

    template<RandomOption Option, uint64_t FixedSeed>
    Array<int> Random<Option, FixedSeed>::random_int(size_t length, int from, int to) {
        Array<int> result(length);
        random_int(result, from, to);
        return result;
    }
    /**
     * Range: [from, to]
     */
    template<RandomOption Option, uint64_t FixedSeed>
    void Random<Option, FixedSeed>::random_int(Array<int>& arr, int from, int to) {
        assert(from <= to && to < INT_MAX);
        if constexpr (HasMKL())
            check_vsl(viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, getInstance(), arr.getLength(), arr.data(), from, to + 1));
        else {
            std::uniform_int_distribution<int> dist(from, to);
            for (int& i : arr)
                i = dist(getInstance());
        }
    }

    template<RandomOption Option, uint64_t FixedSeed>
    auto Random<Option, FixedSeed>::getThreadSeed() const noexcept -> SeedType {
        if constexpr (IsSeedFixed)
            return ThreadPool::isMainThread() ? seed : (seed + ThreadPool::getThreadID() + 1);
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

    template<RandomOption Option, uint64_t FixedSeed>
    constexpr curandRngType_t Random<Option, FixedSeed>::curandRngID() {
    #ifdef PHYSICA_CUDA
        return CURAND_RNG_PSEUDO_MTGP32;
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
