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

#include <forward_list>
#include "clang/Frontend/CompilerInstance.h"
#include "clang/CodeGen/ModuleBuilder.h"
#include "clang/Parse/Parser.h"
#include "llvm/Support/CrashRecoveryContext.h"
#include "ClangImpl/IncrementalAction.h"

namespace Physica::Python {
    /**
     * Reference:
     * [1] clang-repl; https://clang.llvm.org/docs/ClangRepl.html
     */
    class Clang {
        using CompilerInstance = clang::CompilerInstance;
        using CodeGenerator = clang::CodeGenerator;
        using Parser = clang::Parser;
        constexpr static const char* DummyFile = "Unknown";

        struct PartialTranslationUnit {
            clang::TranslationUnitDecl* unitDecl;
            std::unique_ptr<llvm::Module> unitModule;
        };
    private:
        size_t id;
        CompilerInstance ci;
        std::unique_ptr<IncrementalAction> action;
        std::unique_ptr<Parser> parser;
        std::forward_list<PartialTranslationUnit> partialUnitList;
    public:
        Clang();
        Clang(const Clang&) = delete;
        Clang(Clang&&) noexcept = delete;
        ~Clang();
        /* Operators */
        Clang& operator=(const Clang&) = delete;
        Clang& operator=(Clang&&) noexcept = delete;
        /* Operations */
        PartialTranslationUnit& include(const char* path);

        template<class Functor>
        [[nodiscard]] PartialTranslationUnit& makePLU(const char* name, Functor func);
        /* Getters */
        [[nodiscard]] CompilerInstance& getCI() noexcept { return ci; }
        [[nodiscard]] CodeGenerator* getCodeGen() noexcept;
    private:
        /* Operations */
        void makeInvocation();
        void makeOptions();
        void parse();
        void cleanLastUnit() noexcept;
    };

    template<class Functor>
    Clang::PartialTranslationUnit& Clang::makePLU(const char* name, Functor func) {
        using namespace clang;
        Sema& sema = ci.getSema();
        llvm::CrashRecoveryContextCleanupRegistrar<Sema> recoverGuard(&sema);
        Sema::GlobalEagerInstantiationScope GlobalInstantiations(sema, true);
        Sema::LocalEagerInstantiationScope LocalInstantiations(sema);

        PartialTranslationUnit unit{};
        ASTContext& astContext = sema.getASTContext();
        astContext.addTranslationUnitDecl();
        unit.unitDecl = astContext.getTranslationUnitDecl();

        func();

        LocalInstantiations.perform();
        GlobalInstantiations.perform();
        ci.getASTConsumer().HandleTranslationUnit(astContext);

        CodeGenerator* codeGen = getCodeGen();
        if (codeGen != nullptr) {
            auto* unitModule = codeGen->ReleaseModule();
            unitModule->setModuleIdentifier(name);
            unit.unitModule.reset(unitModule);
            codeGen->StartModule(DummyFile, unitModule->getContext());
        }
        partialUnitList.push_front(std::move(unit));
        return partialUnitList.front();
    }
}
