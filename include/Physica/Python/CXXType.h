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

#include <ffi.h>
#include <clang/AST/DeclCXX.h>
#include <nanobind/nanobind.h>

namespace Physica {
    class CXXType {
        clang::CXXRecordDecl* pDecl;
        ffi_type ffiType;

        struct deleter {
            void operator()(void* ptr) { ::operator delete(ptr); }
        };
    public:
        using Ptr = std::unique_ptr<void, deleter>;
    public:
        CXXType() = default;
        CXXType(clang::CXXRecordDecl* pDecl_);
        CXXType(ffi_type ffiType_);
        CXXType(const CXXType&) = default;
        CXXType(CXXType&&) noexcept = default;
        ~CXXType() = default;
        /* Operators */
        CXXType& operator=(CXXType obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] const ffi_type* toFFI() const noexcept { return &ffiType; }
        [[nodiscard]] nanobind::object toPython(void* data) const;
        [[nodiscard]] Ptr allocate() const noexcept;
        void swap(CXXType& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] bool isTrivial() const noexcept { return pDecl == nullptr; }
        [[nodiscard]] size_t getSize() const noexcept { return ffiType.size; }
        [[nodiscard]] size_t getAlign() const noexcept { return ffiType.alignment; }
    };
}
