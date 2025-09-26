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
#include "Physica/Core/Parallel/ThreadPool.h"

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

RandomBase::RandomBase(SeedType seed_, RandomOption option) noexcept : seed(seed_) {
    try {
        reseed(seed_, option);
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

void RandomBase::reseed(SeedType seed_, RandomOption option) {
    seed = seed_;
    SeedType temp = seed;
    if (!ThreadPool::isMainThread())
        for (int i = 0; i < (ThreadPool::getThreadID() + 1) * tseed.size(); ++i)
            RandomSeed::toNextSeed(temp);

    for (auto& elem : tseed) {
        elem = temp;
        RandomSeed::toNextSeed(temp);
    }

    if constexpr (HasMKL())
        check_vsl(vslNewStream(&pStream, rngID_MKL(option), tseed[MKL]));

    if constexpr (HasCUDA()) {
#ifdef PHYSICA_CUDA
        check(curandCreateGenerator(&curand, rngID_cuRAND(option)));
        check(curandSetGeneratorOrdering(curand, CURAND_ORDERING_PSEUDO_DYNAMIC));
        check(curandSetPseudoRandomGeneratorSeed(curand, tseed[CUDA]));
#endif
    }
}
