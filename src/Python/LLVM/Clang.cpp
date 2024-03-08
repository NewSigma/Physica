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
#include <iostream>
#include "clang/Basic/Diagnostic.h"
#include "clang/CodeGen/ObjectFilePCHContainerOperations.h"
#include "clang/Driver/Driver.h"
#include "clang/Driver/Compilation.h"
#include "clang/Driver/Tool.h"
#include "clang/Frontend/TextDiagnosticBuffer.h"
#include "clang/Frontend/TextDiagnosticPrinter.h"
#include "clang/Lex/PreprocessorOptions.h"
#include "llvm/TargetParser/Host.h"
#include "Physica/Python/LLVM/Clang.h"
#include "Physica/Python/Exception/LLVMException.h"
#include "Physica/Macro.h"

namespace Physica::Python {
    Clang::Clang() : ci() {
        using namespace clang;
        {
            auto& PCHOps = *ci.getPCHContainerOperations();
            PCHOps.registerWriter(std::make_unique<ObjectFilePCHContainerWriter>());
            PCHOps.registerReader(std::make_unique<ObjectFilePCHContainerReader>());
        }
        {
            auto* diagOptions = new DiagnosticOptions();
            auto* diagBuffer = new TextDiagnosticPrinter(llvm::outs(), diagOptions);
            auto* diag = new DiagnosticsEngine(new DiagnosticIDs(), diagOptions, diagBuffer, true);
            ci.setDiagnostics(diag);
            makeInvocation();
        }
        auto* memBuffer = llvm::MemoryBuffer::getMemBuffer("").release();
        ci.getPreprocessorOpts().addRemappedFile(DummyFile, memBuffer);
        makeOptions();

        ci.LoadRequestedPlugins();
    }

    void Clang::makeInvocation() {
        const std::string mainExecutableName = llvm::sys::fs::getMainExecutable(nullptr, nullptr);
        std::vector<const char*> args{};
        /* Make args */ {
            args.push_back(mainExecutableName.c_str());
            if constexpr (false) {
                args.push_back("-xcuda");
                args.push_back("--cuda-host-only");
                args.push_back(PHYSICA_CUDA_ARCHITECTURES);
            }
            else
                args.push_back("-xc++");
            args.push_back("-Xclang");
            args.push_back("-fincremental-extensions");
            args.push_back("--compile");
            args.push_back(DummyFile);
            args.push_back(PHYSICA_INCLUDE_DIR);
        }
        using namespace clang;
        auto pDiagOpts = CreateAndPopulateDiagOpts(args);
        auto& diagBuffer = *new TextDiagnosticPrinter(llvm::outs(), pDiagOpts.get());
        DiagnosticsEngine diags(new DiagnosticIDs(), pDiagOpts.release(), &diagBuffer);
        driver::Driver llvmDirver("", llvm::sys::getProcessTriple(), diags);
        std::unique_ptr<driver::Compilation> compilation;
        /* Make compilation */ {
            llvmDirver.setCheckInputsExist(false);
            auto* pCompile = llvmDirver.BuildCompilation(args);
            compilation.reset(pCompile);
            if (pCompile->getArgs().hasArg(driver::options::OPT_v))
                pCompile->getJobs().Print(llvm::errs(), "\n", false);
        }

        const auto& jobs = compilation->getJobs();
        if ((jobs.size() == 0) || !isa<driver::Command>(*jobs.begin()))
            throw LLVMException(llvm::createStringError(llvm::errc::not_supported,"[Error]: Unable to create a driver job"));

        const driver::Command* cmd = cast<driver::Command>(&(*jobs.begin()));
        if (strcmp("clang", cmd->getCreator().getName()) != 0)
            throw LLVMException(llvm::createStringError(llvm::errc::not_supported, "[Error]: Driver initialization failed"));

        const bool success = CompilerInvocation::CreateFromArgs(ci.getInvocation(), cmd->getArguments(), ci.getDiagnostics());
        if (!success)
            throw LLVMException(llvm::createStringError(llvm::errc::not_supported, "[Error]: Failed to create compiler invocation"));
    }

    void Clang::makeOptions() {
        using namespace clang;
        /* Infer the builtin include path if unspecified */ {
            auto& options = ci.getHeaderSearchOpts();
            if (options.UseBuiltinIncludes && options.ResourceDir.empty())
                options.ResourceDir = LLVM_INCLUDE_DIR;
        }
        /* Make target */ {
            auto* target = TargetInfo::CreateTargetInfo(ci.getDiagnostics(), ci.getInvocation().TargetOpts);
            if (target == nullptr)
                throw LLVMException(llvm::createStringError(llvm::errc::not_supported, "[Error]: Target is missing"));
            target->adjust(ci.getDiagnostics(), ci.getLangOpts());
            ci.setTarget(target);
        }
        /* CodeGen */ {
            auto& options = ci.getCodeGenOpts();
            options.ClearASTBeforeBackend = false;
            options.DisableFree = false;
        }
        ci.getFrontendOpts().DisableFree = false;
    }
}
