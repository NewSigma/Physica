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

#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <string_view>

namespace Physica {
    inline std::string readFile(const std::filesystem::path& path) {
        std::ifstream fin(path);
        if (!fin.is_open())
            return {};
        std::string result((std::istreambuf_iterator<char>(fin)), std::istreambuf_iterator<char>());
        while (!result.empty() && std::isspace(result.back()))
            result.pop_back();
        return result;
    }

    inline size_t parseSize(const std::string& text) {
        if (text.empty())
            return 0;
        char* end = nullptr;
        const size_t value = std::strtoull(text.c_str(), &end, 10);
        if (end == text.c_str())
            return 0;
        /* Kernel reports cache sizes with a decimal unit suffix, e.g., "48K" */
        int shift = 0;
        switch (std::toupper(*end)) {
        case 'K':
            shift = 10;
            break;
        case 'M':
            shift = 20;
            break;
        case 'G':
            shift = 30;
            break;
        default:
            break;
        }
        return value << shift;
    }

    inline std::string readCacheAttr(int cpu, int level, std::string_view type, std::string_view attr) {
        const std::filesystem::path root = std::filesystem::path("/sys/devices/system/cpu") / ("cpu" + std::to_string(cpu)) / "cache";
        std::error_code error;
        const auto iter = std::filesystem::directory_iterator(root, error);
        if (error)
            return {};
        /* Directory iteration order is unspecified; always report the lowest cache index */
        constexpr std::string_view Prefix = "index";
        int lowest = std::numeric_limits<int>::max();
        std::string result;
        const std::string levelText = std::to_string(level);
        for (const auto& entry : iter) {
            const std::string name = entry.path().filename().string();
            if (!name.starts_with(Prefix))
                continue;
            if (readFile(entry.path() / "level") != levelText || readFile(entry.path() / "type") != type)
                continue;
            const int index = std::atoi(name.c_str() + Prefix.size());
            if (index < lowest) {
                lowest = index;
                result = readFile(entry.path() / attr);
            }
        }
        return result;
    }

    inline size_t readCacheSize(int cpu, int level, std::string_view type) {
        return parseSize(readCacheAttr(cpu, level, type, "size"));
    }

    inline size_t readCacheLineSize(int cpu) {
        const std::string text = readCacheAttr(cpu, 1, "Data", "coherency_line_size");
        if (text.empty())
            return 0;
        return static_cast<size_t>(std::strtoull(text.c_str(), nullptr, 10));
    }
}
