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

#include <exception>
#include "llvm/Support/Error.h"

namespace Physica::Python {
    class LLVMException : public std::exception {
        std::string msg;
    public:
        explicit LLVMException(const char* err) : msg(err) {}
        explicit LLVMException(llvm::Error err) : msg(llvm::toString(std::move(err))) {}
        ~LLVMException() noexcept override = default;
        const char* what() const noexcept override { return msg.c_str(); }
    };

    inline void llvmCheck(llvm::Error err) {
        if (err) [[unlikely]]
            throw LLVMException(std::move(err));
    }

    template<typename T>
    inline T llvmCheck(llvm::Expected<T>&& E) {
        llvmCheck(E.takeError());
        return T(std::move(*E));
    }

    template<typename T>
    inline T& llvmCheck(llvm::Expected<T&>&& E) {
        llvmCheck(E.takeError());
        return *E;
    }
}
