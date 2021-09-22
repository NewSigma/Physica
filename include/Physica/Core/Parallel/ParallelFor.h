/*
 * Copyright 2021 WeiBo He.
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

#include <vector>
#include <cassert>
#include <cstdlib>
#include "ThreadPool.h"

namespace Physica::Core::Parallel {
    template<class Functor>
    std::vector<std::future<void>> parallel_for(Functor func, size_t loopCount, size_t core) {
        assert(loopCount >= core);
        assert(core > 0);
        const size_t maxLoopPerCore = (loopCount + core - 1) / core;
        size_t from = 0; 
        size_t to = maxLoopPerCore;
        std::vector<std::future<void>> result{};
        result.reserve(core);
        for (size_t i = 0; i < core; ++i) {
            result.push_back(ThreadPool::getInstance().schedule(
                [=](size_t from, size_t to) -> void {
                    for (size_t i = from; i < to; ++i)
                        func(i);
                }, from, to));
            from += maxLoopPerCore;
            const size_t next_to = to + maxLoopPerCore;
            to = next_to > loopCount ? loopCount : next_to;
        }
        return std::move(result);
    }
}
