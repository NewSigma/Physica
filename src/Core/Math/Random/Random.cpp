/*
 * Copyright 2025-2026 Weibo He.
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
#include "Physica/Core/Parallel/MPI.h"
#include "Physica/Core/Utils/Cycler.h"

using namespace Physica;
using namespace Physica::Internal;

RandomState::RandomState(uint64_t seed_, RandomOption option) noexcept : seed(seed_) {
    try {
        reseed(seed_, option);
    }
    catch (std::exception& e) {
        std::cerr << std::format("RNG init failed: {}\n", e.what());
        std::abort();
    }
}

RandomState::~RandomState() {
    if (pStream)
        check_vsl(vslDeleteStream(&pStream));
#ifdef PHYSICA_CUDA
    if (curand)
        check(curandDestroyGenerator(curand));
#endif
}

void RandomState::reseed(uint64_t seed_, RandomOption option) {
    if (seed_ == Dynamic) {
        if constexpr (HostDevAttr::HasRDSEED)
            RandomSeed::rdseed(seed_);
        else if constexpr (HostDevAttr::HasRDRAND)
            RandomSeed::rdrand(seed_);
        else
            seed_ = Cycler::now();
    }

    if constexpr (HasMPI())
        std::ignore = MPI::getInstance();

    seed = seed_;
    seq = SeedSequence<4>({static_cast<uint32_t>(seed),
                           static_cast<uint32_t>(seed >> 32UL),
                           static_cast<uint32_t>(ThreadPool::getThreadID()),
                           static_cast<uint32_t>(MPI::getRank()),});

    if constexpr (HasMKL()) {
        int id_mkl = rngID_MKL(option);
        if (id_mkl != 0) {
            if (pStream) {
                [[maybe_unused]] auto* oldStream = std::exchange(pStream, nullptr);
                check_vsl(vslDeleteStream(&oldStream));
            }
            check_vsl(vslNewStream(&pStream, id_mkl, seq.generate<uint64_t>()));
        }
    }

    if constexpr (HasCUDA()) {
    #ifdef PHYSICA_CUDA
        auto id_curand = rngID_cuRAND(option);
        if (id_curand != CURAND_RNG_TEST) {
            if (curand)
                check(curandDestroyGenerator(std::exchange(curand, nullptr)));
            check(curandCreateGenerator(&curand, id_curand));
            check(curandSetGeneratorOrdering(curand, CURAND_ORDERING_PSEUDO_DYNAMIC));
            check(curandSetPseudoRandomGeneratorSeed(curand, seq.generate<uint64_t>()));
        }
    #endif
    }
}

template<RandomOption Option>
RandomGenerator<Option>::RandomGenerator(uint64_t seed, RandomOption option) noexcept : state(seed, option), gen(state.getSequence()) {}

template<RandomOption Option>
void RandomGenerator<Option>::reseed(uint64_t seed_) {
    state.reseed(seed_, Option);
    gen.seed(state.getSequence());
}
/**
 * Range: [from, to]
 */
template<RandomOption Option>
Array<int> RandomGenerator<Option>::random_int(size_t length, int from, int to) {
    Array<int> result(length);
    random_int(result, from, to);
    return result;
}

template<RandomOption Option>
void RandomGenerator<Option>::random_int(Array<int>& arr, int from, int to) {
    assert(from <= to && to < INT_MAX);
    if constexpr (MKL_Ready)
        check_vsl(viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, state.getMKLStream(), arr.getLength(), arr.data(), from, to + 1));
    else {
        std::uniform_int_distribution<int> dist(from, to);
        for (int& i : arr)
            i = dist(*this);
    }
}

namespace Physica::Internal {
    template class RandomGenerator<PCG32DXSM>;
    template class RandomGenerator<PCG64DXSM>;
    template class RandomGenerator<MCG>;
    template class RandomGenerator<MT19937>;
}
