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
#pragma once

#include <array>
#include <concepts>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <benchmark/benchmark.h>
#include "Utils/Timer.h"
#include "Physica/Core/Math/Algebra/Backend.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/Core/Utils/Builtin.h"
#include "Physica/Core/Utils/Empty.h"
#include "Physica/Core/Utils/MetaProgramming.h"

namespace Physica {
    size_t makeVectorSize(int64_t level, size_t sizeElem) noexcept;
    size_t makeMatrixSize(int64_t level, size_t sizeElem) noexcept;

    namespace Internal {
        template<Backend B, std::invocable F>
        struct BackendEntry {
            constexpr static Backend Kind = B;

            F fn;
        };

        constexpr auto toTuple(auto&& entry) {
            using Entry = decltype(entry);
            if constexpr (instanceof_x<Entry, BackendEntry>)
                return std::make_tuple(std::forward<Entry>(entry));
            else
                return std::tuple<>{};
        }

        template<class Timer, size_t I = 0>
        [[nodiscard]] double time_at(size_t index, Timer& timer, const instanceof<std::tuple> auto& fns) {
            if constexpr (I == std::tuple_size_v<std::remove_cvref_t<decltype(fns)>>)
                unreachable("Backend index is out of range");
            else {
                switch (index) {
                case I:
                    timer.tic();
                    if constexpr (std::same_as<Timer, HostTimer>) {
                        benchmark::DoNotOptimize(std::get<I>(fns).fn());
                        benchmark::ClobberMemory();
                    }
                    else {
                        std::get<I>(fns).fn();
                    }
                    return timer.toc() * 1E9;
                default:
                    return time_at<Timer, I + 1>(index, timer, fns);
                }
            }
        }
    }

    template<Backend B, bool Enabled = true>
    constexpr auto backend(std::invocable auto fn) {
        if constexpr (Enabled)
            return Internal::BackendEntry<B, decltype(fn)>{fn};
        else
            return Empty{};
    }

    void rotate_backends(benchmark::State& state, Either<HostTimer, DeviceTimer> auto&& timer, auto&&... entries) {
        const auto fns = std::tuple_cat(Internal::toTuple(std::forward<decltype(entries)>(entries))...);
        using Tuple = std::remove_cvref_t<decltype(fns)>;
        constexpr size_t NumBackends = std::tuple_size_v<Tuple>;
        static_assert(2 <= NumBackends, "[Error]: rotate_backends needs at least two backends");

        constexpr size_t BaseIdx = [NumBackends]<size_t... I>(std::index_sequence<I...>) consteval noexcept {
            constexpr std::array<Backend, NumBackends> kinds{std::tuple_element_t<I, Tuple>::Kind...};
            for (size_t i = 0; i < NumBackends; ++i)
                if (kinds[i] == Backend::Base)
                    return i;
            return NumBackends;
        }(std::make_index_sequence<NumBackends>{});
        static_assert(BaseIdx != NumBackends, "[Error]: A Base backend is required as reference");

        const auto names = []<size_t... I>(std::index_sequence<I...>) {
            return std::array{std::string(toString(std::tuple_element_t<I, Tuple>::Kind))...};
        }(std::make_index_sequence<NumBackends>{});

        DenseVector<float64, NumBackends> sums{};
        size_t rotation = 0;
        for (auto _ : state) {
            DenseVector<float64, NumBackends> ns{};
            for (size_t i = 0; i < NumBackends; ++i) {
                size_t index = (rotation + i) % NumBackends;
                ns[index] = Internal::time_at(index, timer, fns);
            }

            for (auto&& [sum, elem] : zip(sums, ns))
                sum += elem;
            rotation = (rotation + 1) % NumBackends;

            state.SetIterationTime(ns.sum().toMachine() * 1E-9);
        }
        const double inv = 1.0 / double(state.iterations());
        for (size_t j = 0; j < NumBackends; ++j) {
            state.counters[names[j] + "_ns"] = benchmark::Counter((sums[j] * inv).toMachine());
            if (j != BaseIdx)
                state.counters[names[j] + "_speedup"] = benchmark::Counter((sums[BaseIdx] / sums[j]).toMachine());
        }
    }
}

// Prevent the callee of interest from being inlined, facilitating subsequent analysis at the LLVM IR level.
#ifdef PHYSICA_EmitLLVM
    #define PHYSICA_BENCH(x) [[clang::noinline]] x
#else
    #define PHYSICA_BENCH(x) x
#endif
