/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
    void copy2D(benchmark::State& state) {
        using Buffer = Array<std::byte, Dynamic, HostAllocator<std::byte, 64>>;
        const size_t length = state.range(0);
        const size_t stride = length * state.range(1);
        const size_t span = state.range(2);
        const size_t numRow = span / stride;
        const size_t size = stride * numRow;
        auto src = Buffer::junk(size);
        auto dst = Buffer::junk(size);

        auto base = [&](auto&&...) -> std::byte* {
            for (size_t row = 0; row < numRow; ++row)
                std::memcpy(dst.data() + row * stride, src.data() + row * stride, length);
            return dst.data();
        };
        auto dispatch = [&](auto&&...) -> std::byte* {
            memcpy2D(dst.data(), stride, src.data(), stride, length, numRow);
            return dst.data();
        };
        rotate_backends(state, HostTimer{},
                        backend(base, "base"),
                        backend(dispatch, "dispatch"));
    }
#endif
}

#ifndef PHYSICA_EmitLLVM
BENCHMARK(copy2D)->Name("memcpy2D")->UseManualTime()->ArgsProduct({
    {8, 48, 64, 512, 4096, 10008, 32768, 262144},
    {1, 2},
    {4L * 1024 * 1024, 64L * 1024 * 1024},
});
#endif
