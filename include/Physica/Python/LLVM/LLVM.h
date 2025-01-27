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

#include "llvm/ExecutionEngine/Orc/ThreadSafeModule.h"
#include "llvm/Support/ManagedStatic.h"

namespace Physica::Python {
    class LLVM final {
        using LLVMContext = llvm::LLVMContext;
        using ThreadSafeContext = llvm::orc::ThreadSafeContext;

        llvm::llvm_shutdown_obj shutDownGuard;
        ThreadSafeContext threadSafeContext;
    public:
        LLVM();
        LLVM(const LLVM&) = delete;
        LLVM(LLVM&&) noexcept = delete;
        ~LLVM() = default;
        /* Operators */
        LLVM& operator=(const LLVM&) = delete;
        LLVM& operator=(LLVM&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] ThreadSafeContext& getThreadSafeContext() noexcept { return threadSafeContext; }
        [[nodiscard]] LLVMContext& getContext() noexcept { return *threadSafeContext.getContext(); }
    };
}
