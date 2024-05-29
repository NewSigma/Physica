/*
 * Copyright 2024 WeiBo He.
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

#include "llvm/ExecutionEngine/Orc/LLJIT.h"
#include "clang/Basic/TargetInfo.h"

namespace Physica::Python {
    class Executor {
    public:
        using LLJIT = llvm::orc::LLJIT;
    private:
        std::unique_ptr<LLJIT> jit;
    public:
        Executor() = default;
        Executor(const clang::TargetInfo& targetInfo);
        Executor(const Executor&) = delete;
        Executor(Executor&&) noexcept = default;
        ~Executor() = default;
        /* Operators */
        Executor& operator=(const Executor&) = delete;
        Executor& operator=(Executor&&) noexcept = default;
        /* Operations */
        /* Getters */
        [[nodiscard]] const LLJIT& getJIT() const noexcept { return *jit; }
        [[nodiscard]] LLJIT& getJIT() noexcept { return *jit; }
    };
}
