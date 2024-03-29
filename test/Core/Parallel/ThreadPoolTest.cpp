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
#include <iostream>
#include "Physica/Core/Parallel/ParallelFor.h"

using namespace Physica::Core::Parallel;

void func([[maybe_unused]] size_t i) {
    printf("Thread ID: %ld\n", ThreadPool::getThreadInfo().id);
}

int main() {
    ThreadPool::initThreadPool(4);
    ThreadPool& pool = ThreadPool::getInstance();
    parallel_for(func, 20, 4);
    pool.shouldExit();
    ThreadPool::deInitThreadPool();
    return 0;
}
