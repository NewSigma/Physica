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
#include "clang/Frontend/CompilerInstance.h"
#include "clang/FrontendTool/Utils.h"
#include "clang/CodeGen/CodeGenAction.h"
#include "clang/Sema/CodeCompleteConsumer.h"
#include "Physica/Python/LLVM/ClangImpl/IncrementalAction.h"
#include "Physica/Python/Exception/LLVMException.h"

namespace Physica::Python {
    IncrementalAction::IncrementalAction(CompilerInstance& ci, llvm::LLVMContext& context)
            : Base(makeAction(ci, context)), isTerminating(false) {}

    void IncrementalAction::ExecuteAction() {
        CompilerInstance& ci = getCompilerInstance();
        assert(ci.hasPreprocessor() && "[Error]: Invalid CI");
        if (hasCodeCompletionSupport() && !ci.getFrontendOpts().CodeCompletionAt.FileName.empty())
            ci.createCodeCompletionConsumer();

        ci.getPreprocessor().EnterMainSourceFile();
        if (!ci.hasSema()) {
            clang::CodeCompleteConsumer* completeConsumer = nullptr;
            if (ci.hasCodeCompletionConsumer())
                completeConsumer = &ci.getCodeCompletionConsumer();
            ci.createSema(getTranslationUnitKind(), completeConsumer);
        }
    }

    void IncrementalAction::EndSourceFile() {
        if (isTerminating && (getWrapped() != nullptr))
            Base::EndSourceFile();
    }

    void IncrementalAction::finalizeAction() {
        assert(!isTerminating && "[Error]: Double finalized");
        if (isTerminating) [[unlikely]]
            return;
        isTerminating = true;
        EndSourceFile();
    }

    std::unique_ptr<clang::FrontendAction> IncrementalAction::makeAction(CompilerInstance& ci, llvm::LLVMContext& context) {
        using namespace clang::frontend;
        std::unique_ptr<FrontendAction> result{};
        switch (ci.getFrontendOpts().ProgramAction) {
        case ASTDump:
            [[fallthrough]];
        case ASTPrint:
            [[fallthrough]];
        case ParseSyntaxOnly:
            result = clang::CreateFrontendAction(ci);
            break;
        case PluginAction:
            [[fallthrough]];
        case EmitAssembly:
            [[fallthrough]];
        case EmitBC:
            [[fallthrough]];
        case EmitObj:
            [[fallthrough]];
        case PrintPreprocessedInput:
            [[fallthrough]];
        case EmitLLVMOnly:
            result.reset(new clang::EmitLLVMOnlyAction(&context));
            break;
        default:
            throw LLVMException(llvm::createStringError(std::errc::state_not_recoverable,
                "[Error]: Incremental mode for action %d is not supported", ci.getFrontendOpts().ProgramAction));
        }
        return result;
    }
}
