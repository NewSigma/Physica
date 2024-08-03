/*
 * Copyright 2024 Weibo He.
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
#include "Physica/Python/LLVM/Clang.h"
#include "Physica/Python/LLVM/Executor.h"
#include "clang/AST/DeclCXX.h"

namespace Physica::Python {
    class PhysicaPython final {
        using This = PhysicaPython;
        using ExecutorAddr = llvm::orc::ExecutorAddr;
        using ClassDeclMap = std::unordered_map<const char*, clang::CXXRecordDecl*>;
        using LLJIT = typename Executor::LLJIT;

        Clang clang;
        Executor exec;
        ClassDeclMap classDeclMap;
    public:
        ~PhysicaPython() = default;
        /* Getters */
        [[nodiscard]] const Clang& getClang() const noexcept { return clang; }
        [[nodiscard]] Clang& getClang() noexcept { return clang; }
        [[nodiscard]] const Executor& getExec() const noexcept { return exec; }
        [[nodiscard]] Executor& getExec() noexcept { return exec; }
        [[nodiscard]] const LLJIT& getJIT() const noexcept { return exec.getJIT(); }
        [[nodiscard]] LLJIT& getJIT() noexcept { return exec.getJIT(); }
        [[nodiscard]] const ClassDeclMap& getClassDeclMap() const noexcept { return classDeclMap; }
        /* Static members */
        [[nodiscard]] static This& getInstance() noexcept;
    private:
        PhysicaPython();
        PhysicaPython(const This&) = default;
        PhysicaPython(This&&) noexcept = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
    };
}
