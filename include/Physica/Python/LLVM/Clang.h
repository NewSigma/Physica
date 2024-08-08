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

#include <forward_list>
#include <llvm/IR/Module.h>
#include "llvm/Support/CrashRecoveryContext.h"
#include "clang/Frontend/CompilerInstance.h"
#include "clang/CodeGen/ModuleBuilder.h"
#include "clang/Parse/Parser.h"
#include "ClangImpl/IncrementalAction.h"

namespace Physica::Python {
    class HeaderManager;
    /**
     * Reference:
     * [1] clang-repl; https://clang.llvm.org/docs/ClangRepl.html
     */
    class Clang {
        using This = Clang;
        using CompilerInstance = clang::CompilerInstance;
        using CodeGenerator = clang::CodeGenerator;
        using Parser = clang::Parser;
        using TranslationUnitDecl = clang::TranslationUnitDecl;
        constexpr static const char* DummyFile = "Unknown";
    public:
        struct PartialTranslationUnit {
            TranslationUnitDecl* unitDecl;
            std::unique_ptr<llvm::Module> unitModule;
        };
    private:
        CompilerInstance ci;
        std::unique_ptr<IncrementalAction> action;
        std::unique_ptr<Parser> parser;
        std::forward_list<PartialTranslationUnit> partialUnitList;
        HeaderManager* pHeaderManager;
    public:
        Clang();
        Clang(const Clang&) = delete;
        Clang(Clang&&) noexcept = delete;
        ~Clang() = default;
        /* Operators */
        Clang& operator=(const Clang&) = delete;
        Clang& operator=(Clang&&) noexcept = delete;
        /* Operations */
        const clang::NamedDecl* include(const char* includeName);
        PartialTranslationUnit& compile(const char* moduleName);
        /* Getters */
        [[nodiscard]] CompilerInstance& getCI() noexcept { return ci; }
        [[nodiscard]] const clang::ASTContext& getASTContext() const noexcept { return ci.getASTContext(); }
        [[nodiscard]] clang::ASTContext& getASTContext() noexcept { return ci.getASTContext(); }
        [[nodiscard]] const CodeGenerator& getCodeGen() const noexcept;
        [[nodiscard]] CodeGenerator& getCodeGen() noexcept;
        [[nodiscard]] Parser& getParser() noexcept { return *parser; }
        [[nodiscard]] const HeaderManager& getHeaderManager() const noexcept { return *pHeaderManager; }
    private:
        /* Operations */
        void makeInvocation();
        void makeOptions();
        void parse();
        void cleanLastUnit() noexcept;
    };
}
