/*
 * Copyright 2023-2026 Weibo He.
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
#include "Physica/Core/Utils/Container/Array.h"
#include "Generator/PCG.h"
#include "SeedSequence.h"

namespace Physica {
    enum RandomOption : char {
        PCG32DXSM,
        PCG64DXSM,
        MCG,
        MT19937,
    };

    namespace Internal {
        class PHYSICA_API RandomBase {
            struct Private; // Do not conflict with other pointers
        protected:
        #ifdef PHYSICA_MKL
            using mkl_handle = VSLStreamStatePtr;
        #else
            using mkl_handle = void*;
        #endif

        #ifdef PHYSICA_CUDA
            using curand_handle = curandGenerator_t;
        #else
            using curand_handle = Private*;
        #endif

            SeedSequence<> seq;
            mkl_handle pStream = nullptr;
            [[no_unique_address]] curand_handle curand = nullptr;
        private:
            uint64_t seed{};
        public:
            RandomBase() = delete;
            ~RandomBase();
            /* Getters */
            [[nodiscard]] uint64_t getSeed() const noexcept { return seed; }
        protected:
            RandomBase(uint64_t seed_, RandomOption option) noexcept;
            /* Operations */
            void reseed(RandomOption option);
            void reseed(uint64_t seed_, RandomOption option);
        };
    }

    constexpr int rngID_MKL(RandomOption option) noexcept {
    #ifdef PHYSICA_MKL
        switch (option) {
        case MCG:
            return VSL_BRNG_MCG31;
        case MT19937:
            return VSL_BRNG_MT19937;
        default:
            return 0;
        }
    #else
        return 0;
    #endif
    }

    constexpr auto rngID_cuRAND([[maybe_unused]] RandomOption option) noexcept {
    #ifdef PHYSICA_CUDA
        switch (option) {
        case MT19937:
            return CURAND_RNG_PSEUDO_MTGP32;
        default:
            static_assert(CURAND_RNG_TEST == 0, "Returns 0 by default");
            return CURAND_RNG_TEST;
        }
    #else
        return 0;
    #endif
    }
    /**
     * \class Random provides a general, per-thread, reusable random generator implementation.
     * We use the same random number generator as NumPy [1].
     *
     * Note: Even if we fix the random seed, thread stealing may still break reproducibility.
     *
     * \tparam FixedSeed:
     * Multi-threading will break reproducibility and fixing random seed is not enough.
     * Therefore, \tparam FixedSeed is declared as a template param,
     * making it possible for other parts of the program to check at compile time whether the seed is fixed.
     *
     * Reference:
     * [1] NumPy BitGenerators; https://numpy.org/doc/stable/reference/random/bit_generators/index.html
     */
    template<RandomOption Option = PCG64DXSM, uint64_t FixedSeed = Physica::Dynamic>
    class PHYSICA_API Random : public Internal::RandomBase {
        using This = Random<Option, FixedSeed>;
        using Base = Internal::RandomBase;
        using typename Base::mkl_handle;
        using typename Base::curand_handle;

        template<RandomOption> struct GeneratorImpl;
        template<> struct GeneratorImpl<PCG32DXSM> { using Type = Internal::setseq_base<uint32_t, uint64_t, Internal::DXSM>; };
        template<> struct GeneratorImpl<PCG64DXSM> { using Type = Internal::setseq_base<uint64_t, __uint128_t, Internal::DXSM>; };
        template<> struct GeneratorImpl<MCG> { using Type = std::minstd_rand; };
        template<> struct GeneratorImpl<MT19937> { using Type = std::mt19937; };

        using Generator = GeneratorImpl<Option>::Type;
    public:
        using result_type = Generator::result_type;
        constexpr static bool IsSeedFixed = Traits<This>::IsSeedFixed;
        constexpr static bool MKL_Ready = rngID_MKL(Option) != 0;
        constexpr static bool cuRAND_Ready = (int)rngID_cuRAND(Option) != 0;
    private:
        Generator gen;
    public:
        ~Random() = default;
        /* Operators */
        [[nodiscard]] result_type operator()() { return gen(); }
        [[nodiscard]] operator mkl_handle() noexcept requires(MKL_Ready);
        [[nodiscard]] operator curand_handle() noexcept requires(cuRAND_Ready);
        /* Operations */
        void reseed();
        void reseed(uint64_t seed_) requires(!IsSeedFixed);
        /* Getters */
        [[nodiscard]] Generator& generator() noexcept { return gen; }
        /* Static members */
        [[nodiscard]] constexpr static result_type min() { return Generator::min(); }
        [[nodiscard]] constexpr static result_type max() { return Generator::max(); }

        [[nodiscard]] static This& getInstance() noexcept;
        [[nodiscard]] static bool coin() noexcept;
        [[nodiscard]] static int select(int size) noexcept;
        [[nodiscard]] static Array<int> random_int(size_t length, int from, int to);
        static void random_int(Array<int>& arr, int from, int to);
    private:
        Random() noexcept;
        Random(const This&) = default;
        Random(This&&) noexcept = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
    };

    template<RandomOption Option, uint64_t FixedSeed>
    Random<Option, FixedSeed>::Random() noexcept : Base(FixedSeed, Option) {
        gen.seed(seq);
    }

    template<RandomOption Option, uint64_t FixedSeed>
    Random<Option, FixedSeed>::operator mkl_handle() noexcept requires(MKL_Ready){
        static_assert(MKL_Ready, "[Error]: MKL does not support this RNG");
        return Base::pStream;
    }

    template<RandomOption Option, uint64_t FixedSeed>
    Random<Option, FixedSeed>::operator curand_handle() noexcept requires(cuRAND_Ready) {
        static_assert(cuRAND_Ready, "[Error]: cuRAND does not support this RNG");
        return Base::curand;
    }

    template<RandomOption Option, uint64_t FixedSeed>
    void Random<Option, FixedSeed>::reseed() {
        if constexpr (IsSeedFixed)
            Base::reseed(FixedSeed, Option);
        else
            Base::reseed(Option);
        gen.seed(seq);
    }

    template<RandomOption Option, uint64_t FixedSeed>
    void Random<Option, FixedSeed>::reseed(uint64_t seed) requires(!IsSeedFixed) {
        Base::reseed(seed, Option);
        gen.seed(seq);
    }

    template<RandomOption Option, uint64_t FixedSeed>
    auto Random<Option, FixedSeed>::getInstance() noexcept -> This& {
        thread_local static This instance{};
        return instance;
    }

    template<RandomOption Option, uint64_t FixedSeed>
    bool Random<Option, FixedSeed>::coin() noexcept {
        return select(2);
    }

    template<RandomOption Option, uint64_t FixedSeed>
    int Random<Option, FixedSeed>::select(int size) noexcept {
        assert(size > 0 && "[Error]: select size must be positive");
        return std::uniform_int_distribution<>(0, size - 1)(getInstance());
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
        if constexpr (MKL_Ready)
            check_vsl(viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, getInstance(), arr.getLength(), arr.data(), from, to + 1));
        else {
            std::uniform_int_distribution<int> dist(from, to);
            for (int& i : arr)
                i = dist(getInstance());
        }
    }

    template<class T>
    concept RNG = instanceof_xx<T, Random>;
}

namespace Physica {
    template<RandomOption Option, uint64_t FixedSeed>
    class Traits<Random<Option, FixedSeed>> {
    public:
        constexpr static bool IsSeedFixed = FixedSeed != Dynamic;
    };
}
