/*
 * Copyright 2020-2025 Weibo He.
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
#include "LoggerType.h"

namespace Physica {
    template<size_t N>
    constexpr size_t getArgCount(const char (&format)[N]);

    template<size_t N>
    constexpr ArgType getArgType(const char (&format)[N], int index);

    template<size_t ArgCount, size_t N>
    constexpr std::array<ArgType, ArgCount> analyzeFormatString(const char(&format)[N]);

    template<size_t N, std::size_t... Indices>
    constexpr std::array<ArgType, sizeof...(Indices)>
    analyzeFormatStringHelper(const char (&format)[N], std::index_sequence<Indices...>);
}

#include "FormatAnalyzerImpl.h"
