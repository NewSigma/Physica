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
#include <format>
#include "Physica/Python/LLVM/ASTCursor.h"

using namespace Physica;

ASTCursor::ASTCursor(clang::ASTContext& ctx) {
    push(ctx.getTranslationUnitDecl());
}

std::string ASTCursor::toString() const {
    return formatDecl(astStack.top());
}

void ASTCursor::push(const DeclContext* decl) {
    astStack.push(decl);
    traverse();
}

void ASTCursor::pop() {
    if (astStack.size() == 1)
        return;
    astStack.pop();
    traverse();
}

void ASTCursor::ls() const {
    for (size_t i = 0; i < childs.getLength(); ++i) {
        llvm::outs() << i << '\t';

        const auto* decl = childs[i];
        if (llvm::isa<DeclContext>(decl)) {
            const auto* p = llvm::cast<DeclContext>(decl);
            llvm::outs() << std::format("+ {}\n", formatDecl(p));
        }
        else
            decl->dump();
    }
}

std::string ASTCursor::cd(size_t index) {
    if (!(index < size()))
        throw std::out_of_range("[Error]: Invalid child");

    const auto* decl = childs[index];
    if (!llvm::isa<DeclContext>(decl))
        throw std::invalid_argument("[Error]: This is not a DeclContext");
    push(llvm::cast<DeclContext>(decl));
    return toString();
}

void ASTCursor::dump() const {
    astStack.top()->dumpAsDecl();
}

void ASTCursor::reset() {
    while (astStack.size() > 1)
        astStack.pop();
    traverse();
}

void ASTCursor::swap(This& __restrict obj) noexcept {
    assert(this != &obj && "[Error]: Self swap is likely a bug");
    astStack.swap(obj.astStack);
    childs.swap(obj.childs);
}

void ASTCursor::traverse() {
    childs.clear();
    for (const auto* decl : astStack.top()->decls())
        childs.append(decl);
}

std::string ASTCursor::formatDecl(const DeclContext* decls) {
    if (llvm::isa<clang::NamedDecl>(decls)) {
        const auto* p = llvm::cast<clang::NamedDecl>(decls);
        return std::format("{}({})", p->getNameAsString(), p->getDeclKindName());
    }
    else
        return std::format("{}", decls->getDeclKindName());
}
