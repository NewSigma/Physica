/*
 * Copyright 2026 Weibo He.
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

#include <exception>
#include <ranges>
#include "Physica/Core/Parallel/Task/Task.h"

namespace Physica {
    [[nodiscard]] Task when_all(std::ranges::range auto&& tasks) {
        std::exception_ptr ex = nullptr;
        for (auto& task : decay_rvalue(std::forward<decltype(tasks)>(tasks))) {
            try {
                co_await task;
            }
            catch (...) {
                if (ex == nullptr)
                    ex = std::current_exception();
            }
        }

        if (ex) [[unlikely]]
            std::rethrow_exception(ex);
    }
}
