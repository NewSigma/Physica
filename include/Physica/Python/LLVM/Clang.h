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

#include "clang/Frontend/CompilerInstance.h"

namespace Physica::Python {
    class Clang {
        using CompilerInstance = clang::CompilerInstance;
        constexpr static const char* DummyFile = "<<< inputs >>>";

        CompilerInstance ci;
    public:
        Clang();
        Clang(const Clang&) = delete;
        Clang(Clang&&) noexcept = delete;
        ~Clang() = default;
        /* Operators */
        Clang& operator=(const Clang&) = delete;
        Clang& operator=(Clang&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] CompilerInstance& getCI() noexcept { return ci; }
    private:
        /* Operations */
        void makeInvocation();
        void makeOptions();
    };
}
