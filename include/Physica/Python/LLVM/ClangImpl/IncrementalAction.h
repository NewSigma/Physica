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

#include "llvm/IR/LLVMContext.h"
#include "clang/Frontend/FrontendAction.h"

namespace llvm {
    class GlobalValue;
}

namespace clang::CodeGen {
    /**
     * Hack clang CodeGen
     */
    class CodeGenModule {
    public:
        void EmitGlobalDefinition(GlobalDecl D, llvm::GlobalValue *GV = nullptr);
    };
}

namespace Physica {
    class IncrementalAction : public clang::WrapperFrontendAction {
        using This = IncrementalAction;
        using Base = clang::WrapperFrontendAction;
        using CompilerInstance = clang::CompilerInstance;
        using LLVMContext = llvm::LLVMContext;
    public:
        IncrementalAction(CompilerInstance& ci, llvm::LLVMContext& context);
        IncrementalAction(const This&) = delete;
        IncrementalAction(This&&) = delete;
        ~IncrementalAction() override;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void ExecuteAction() override;
        void EndSourceFile() override;
        /* Getters */
        [[nodiscard]] clang::FrontendAction* getWrapped() const { return WrappedAction.get(); }
        [[nodiscard]] clang::TranslationUnitKind getTranslationUnitKind() override { return clang::TU_Incremental; }
    private:
        static std::unique_ptr<clang::FrontendAction> makeAction(CompilerInstance& ci, llvm::LLVMContext& context);
    };
}
