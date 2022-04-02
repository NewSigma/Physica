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
    std::vector<std::future<void>> parallel_for(Functor func, unsigned int loopCount, unsigned int core) {
        assert(loopCount >= core);
        assert(core > 0);
        const unsigned int maxLoopPerCore = (loopCount + core - 1) / core;
        unsigned int from = 0; 
        unsigned int to = maxLoopPerCore;
        std::vector<std::future<void>> result{};
        result.reserve(core);
        for (unsigned int i = 0; i < core; ++i) {
            result.push_back(ThreadPool::getInstance().schedule(
                [=](unsigned int from_, unsigned int to_) -> void {
                    for (unsigned int i = from_; i < to_; ++i)
                        func(i);
                }, from, to));
            from += maxLoopPerCore;
            const unsigned int next_to = to + maxLoopPerCore;
            to = next_to > loopCount ? loopCount : next_to;
        }
        return result;
    }
}
