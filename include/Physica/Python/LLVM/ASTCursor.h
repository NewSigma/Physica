/*
 * Copyright 2025 Weibo He.
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

#include <stack>
#include "clang/AST/ASTContext.h"
#include "Physica/Core/Utils/Container/Array.h"

namespace Physica::Python {
    class ASTCursor {
        using This = ASTCursor;
        using DeclContext = clang::DeclContext;
        using Decl = clang::Decl;

        std::stack<const DeclContext*> astStack;
        Array<const Decl*> childs;
    public:
        ASTCursor(clang::ASTContext& ctx_);
        ASTCursor(const This&) = default;
        ASTCursor(This&&) noexcept = default;
        ~ASTCursor() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] std::string toString() const;

        void push(const DeclContext* decl);
        void pop();

        void ls() const;
        [[nodiscard]] std::string cd(size_t index);
        void dump() const;
        void reset();
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t size() const noexcept { return childs.getLength(); }
    private:
        void traverse();
        static std::string formatDecl(const DeclContext* decls);
    };
}
