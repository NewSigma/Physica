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

#include <llvm/Support/Error.h>

namespace Physica {
    class LLVMException : public std::runtime_error {
        using Base = std::runtime_error;
    public:
        explicit LLVMException(const char* err) : Base(err) {}
        explicit LLVMException(llvm::Error err) : Base(llvm::toString(std::move(err))) {}
    };

    inline void check(llvm::Error err) {
        if (err) [[unlikely]]
            throw LLVMException(std::move(err));
    }

    template<typename T>
    T check(llvm::Expected<T>&& E) {
        check(E.takeError());
        return T(std::move(*E));
    }

    template<typename T>
    T& check(llvm::Expected<T&>&& E) {
        check(E.takeError());
        return *E;
    }
}
