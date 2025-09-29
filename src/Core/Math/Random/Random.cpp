/*
 * Copyright 2025 Weibo He.
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
#include <iostream>
#include <format>
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Math/Random/RandomSeed.h"
#include "Physica/Core/Parallel/ThreadPool.h"
#ifdef PHYSICA_MPI
    #include "Physica/Core/Parallel/MPIContext.h"
#endif

using namespace Physica;
using namespace Physica::Internal;

namespace {
    [[maybe_unused]] int rngID_MKL(RandomOption) noexcept {
    #ifdef PHYSICA_MKL
        return VSL_BRNG_MT19937;
    #else
        return 0;
    #endif
    }

    [[maybe_unused]] curandRngType_t rngID_cuRAND(RandomOption) noexcept {
    #ifdef PHYSICA_CUDA
        return CURAND_RNG_PSEUDO_MTGP32;
    #else
        return 0;
    #endif
    }
}

RandomBase::RandomBase(RandomOption option) noexcept {
    try {
        reseed(option);
    }
    catch (std::exception& e) {
        std::cerr << std::format("RNG init failed: {}\n", e.what());
        std::abort();
    }
}

RandomBase::~RandomBase() {
    check_vsl(vslDeleteStream(&pStream));
#ifdef PHYSICA_CUDA
    check(curandDestroyGenerator(curand));
#endif
}

void RandomBase::reseed(RandomOption option) {
    SeedType seed{};
    RandomSeed::rdseed(seed);
    reseed(seed, option);
}

void RandomBase::reseed(SeedType seed_, RandomOption option) {
    seed = seed_;
    seq = SeedSequence<>({
        static_cast<uint32_t>(seed),
        static_cast<uint32_t>(seed >> 32UL),
        static_cast<uint32_t>(ThreadPool::getThreadID()),
    #ifdef PHYSICA_MPI
        static_cast<uint32_t>(HasMPI() ? MPIContext::getInstance().getProcessID() : 0)
    #endif
    });

    if constexpr (HasMKL())
        check_vsl(vslNewStream(&pStream, rngID_MKL(option), seq.generate<uint64_t>()));

    if constexpr (HasCUDA()) {
    #ifdef PHYSICA_CUDA
        check(curandCreateGenerator(&curand, rngID_cuRAND(option)));
        check(curandSetGeneratorOrdering(curand, CURAND_ORDERING_PSEUDO_DYNAMIC));
        check(curandSetPseudoRandomGeneratorSeed(curand, seq.generate<uint64_t>()));
    #endif
    }
}
