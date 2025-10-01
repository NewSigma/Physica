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
#include "Physica/Core/Parallel/MPIContext.h"

using namespace Physica;
using namespace Physica::Internal;

RandomBase::RandomBase(uint64_t seed_, RandomOption option) noexcept {
    try {
        if (seed_ == Dynamic)
            reseed(option);
        else
            reseed(seed_, option);
    }
    catch (std::exception& e) {
        std::cerr << std::format("RNG init failed: {}\n", e.what());
        std::abort();
    }
}

RandomBase::~RandomBase() {
    if (pStream)
        check_vsl(vslDeleteStream(&pStream));
#ifdef PHYSICA_CUDA
    if (curand)
        check(curandDestroyGenerator(curand));
#endif
}

void RandomBase::reseed(RandomOption option) {
    uint64_t seed{};
    RandomSeed::rdseed(seed);
    reseed(seed, option);
}

void RandomBase::reseed(uint64_t seed_, RandomOption option) {
    seed = seed_;
    seq = SeedSequence<4>({
        static_cast<uint32_t>(seed),
        static_cast<uint32_t>(seed >> 32UL),
        static_cast<uint32_t>(ThreadPool::getThreadID()),
        static_cast<uint32_t>(MPIContext::getInstance().getProcessID())
    });

    
    if constexpr (HasMKL()) {
        int id_mkl = rngID_MKL(option);
        if (id_mkl != 0)
            check_vsl(vslNewStream(&pStream, id_mkl, seq.generate<uint64_t>()));
    }

    if constexpr (HasCUDA()) {
    #ifdef PHYSICA_CUDA
        auto id_curand = rngID_cuRAND(option);
        if (id_curand != CURAND_RNG_TEST) {
            check(curandCreateGenerator(&curand, id_curand));
            check(curandSetGeneratorOrdering(curand, CURAND_ORDERING_PSEUDO_DYNAMIC));
            check(curandSetPseudoRandomGeneratorSeed(curand, seq.generate<uint64_t>()));
        }
    #endif
    }
}
