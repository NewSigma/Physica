/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <cstddef>
#include <cstring>
#include "Physica/Core/Utils/Allocator/HostAllocator.h"
#include "Physica/Core/Utils/Container/Array.h"
#include "Physica/Core/Utils/Memory.h"
#include "Benchmark.h"

using namespace Physica;

namespace {
#ifndef PHYSICA_EmitLLVM
    void set2D(benchmark::State& state) {
        using Buffer = Array<std::byte, Dynamic, HostAllocator<std::byte, 64>>;
        const size_t length = state.range(0);
        const size_t stride = length * state.range(1);
        const size_t span = state.range(2);
        const size_t numRow = span / stride;
        const size_t size = stride * numRow;
        const int value = int(state.range(3));
        auto dst = Buffer::junk(size);

        auto base = [&](auto&&...) -> std::byte* {
            for (size_t row = 0; row < numRow; ++row)
                std::memset(dst.data() + row * stride, value, length);
            return dst.data();
        };
        auto dispatch = [&](auto&&...) -> std::byte* {
            PHYSICA_BENCH(memset2D(dst.data(), stride, value, length, numRow));
            return dst.data();
        };
        rotate_backends(state, HostTimer{},
                        backend(base, "base"),
                        backend(dispatch, "dispatch"));
    }
#endif
}

#ifndef PHYSICA_EmitLLVM
BENCHMARK(set2D)->Name("memset2D")->UseManualTime()->ArgsProduct({
    {8, 48, 64, 512, 4096, 10008, 32768, 262144},
    {1, 2},
    {4L * 1024 * 1024, 64L * 1024 * 1024},
    {0},
});
#endif
