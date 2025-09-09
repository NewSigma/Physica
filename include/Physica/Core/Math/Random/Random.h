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
#include "Physica/Core/Utils/Container/Array.h"
#ifdef PHYSICA_CUDA
    #include "Physica/Core/Exception/CUDA/cuRAND.cuh"
#endif
#include "Physica/Core/Exception/MKL/VSL.h"
#include "RandomSeed.h"

namespace Physica {
    enum RandomOption : char {
        MT19937
    };

    namespace Internal {
        class PHYSICA_API RandomBase {
        public:
            using SeedType = uint64_t;
        private:
            VSLStreamStatePtr pStream = nullptr;
            [[no_unique_address]] curandGenerator_t curand{};

            SeedType seed{};
            SeedType tseed{};
        public:
            RandomBase() = delete;
            ~RandomBase();
            /* Operators */
            [[nodiscard]] operator VSLStreamStatePtr() noexcept { return pStream; }
            [[nodiscard]] operator curandGenerator_t() noexcept { return curand; }
            /* Getters */
            [[nodiscard]] SeedType getSeed() const noexcept { return seed; }
            [[nodiscard]] SeedType getTSeed() const noexcept { return tseed; }
        protected:
            RandomBase(SeedType seed_, RandomOption option, bool fixed) noexcept;
            /* Operations */
            void reseed(SeedType seed_, RandomOption option, bool fixed);
        };

        class QRandomBase {};
    }
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
        using Base = Internal::RandomBase;
    public:
        using result_type = uint64_t;
        constexpr static bool IsSeedFixed = Traits<This>::IsSeedFixed;
    private:
        using GenType = std::mt19937;

        GenType gen;
    public:
        ~Random() = default;
        /* Operators */
        [[nodiscard]] result_type operator()() { return gen(); }
        /* Operations */
        void reseed();
        void reseed(SeedType seed_);
        /* Getters */
        [[nodiscard]] GenType& getGen() noexcept { return gen; }
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
        /* Static members */
        static SeedType genSeed() noexcept;
    };

    template<RandomOption Option, uint64_t FixedSeed>
    Random<Option, FixedSeed>::Random() noexcept : Base(genSeed(), Option, FixedSeed) {
        gen.seed(getTSeed());
    }

    template<RandomOption Option, uint64_t FixedSeed>
    void Random<Option, FixedSeed>::reseed() {
        reseed(genSeed());
    }

    template<RandomOption Option, uint64_t FixedSeed>
    void Random<Option, FixedSeed>::reseed(SeedType seed) {
        if constexpr (IsSeedFixed)
            assert(seed == FixedSeed);
        Base::reseed(seed, Option, FixedSeed);
        gen.seed(getTSeed());
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
    auto Random<Option, FixedSeed>::genSeed() noexcept -> SeedType {
        SeedType seed{};
        if constexpr (IsSeedFixed)
            seed = FixedSeed;
        else
            RandomSeed::rdrand(seed);
        return seed;
    }
}

namespace Physica {
    template<RandomOption Option, uint64_t FixedSeed>
    class Traits<Random<Option, FixedSeed>> {
    public:
        constexpr static bool IsSeedFixed = FixedSeed != Dynamic;
    };
}
