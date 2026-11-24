/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
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
