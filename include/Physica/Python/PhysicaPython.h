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

#include "Physica/Python/LLVM/Clang.h"
#include "Physica/Python/LLVM/Executor.h"

namespace Physica::Python {
    class PhysicaPython final {
        using This = PhysicaPython;

        Clang clang;
        Executor exec;
    public:
        ~PhysicaPython() = default;
        /* Getters */
        [[nodiscard]] Clang& getClang() noexcept { return clang; }
        [[nodiscard]] Executor& getExec() noexcept { return exec; }
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
