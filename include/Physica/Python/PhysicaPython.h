/*
 * Copyright 2024-2025 Weibo He.
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

#include <unordered_map>
#include <filesystem>
#include "llvm/ExecutionEngine/Orc/LLJIT.h"
#include "Physica/Python/LLVM/Clang.h"
#include "CXXType.h"

namespace Physica::Python {
    class PhysicaPython final {
        using This = PhysicaPython;
        using ExecutorAddr = llvm::orc::ExecutorAddr;
        using StrTypeMap = std::unordered_map<std::string, CXXType>;
        using LLJIT = llvm::orc::LLJIT;

        Clang clang;
        std::unique_ptr<LLJIT> jit;
        StrTypeMap strTypeMap;
    public:
        PhysicaPython(std::filesystem::path root_);
        PhysicaPython(const This&) = delete;
        PhysicaPython(This&&) noexcept = delete;
        ~PhysicaPython() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] const CXXType& toCXXType(const std::string& typeName) const { return strTypeMap.at(typeName); }
        [[nodiscard]] const CXXType& toCXXType(py::handle handle) const { return toCXXType(std::string(py::str(handle))); }

        void compile(const char* moduleName);
        /* Getters */
        [[nodiscard]] Clang& getClang() noexcept { return clang; }
        [[nodiscard]] LLJIT& getJIT() noexcept { return *jit; }
        [[nodiscard]] const StrTypeMap& getStrTypeMap() const noexcept { return strTypeMap; }
        /* Static members */
        [[nodiscard]] static This& getInstance() noexcept;
    };
}
