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
#include <sstream>
#include "clang/AST/DeclContextInternals.h"
#include "clang/Basic/Diagnostic.h"
#include "clang/CodeGen/ObjectFilePCHContainerOperations.h"
#include "clang/CodeGen/CodeGenAction.h"
#include "clang/Driver/Driver.h"
#include "clang/Driver/Compilation.h"
#include "clang/Driver/Tool.h"
#include "clang/Frontend/TextDiagnosticBuffer.h"
#include "clang/Frontend/TextDiagnosticPrinter.h"
#include "clang/Lex/PreprocessorOptions.h"
#include "clang/Lex/Preprocessor.h"
#include "clang/Sema/Sema.h"
#include "llvm/TargetParser/Host.h"
#include "llvm/IR/Module.h"
#include "llvm/Support/CrashRecoveryContext.h"
#include "Physica/Python/LLVM/LLVM.h"
#include "Physica/Python/LLVM/Clang.h"
#include "Physica/Python/Exception/LLVMException.h"
#include "Physica/Macro.h"

namespace Physica::Python {
    Clang::Clang() : id(0), ci() {
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

        action = std::make_unique<IncrementalAction>(ci, LLVM::getInstance().getContext());
        ci.ExecuteAction(*action);

        parser = std::make_unique<Parser>(ci.getPreprocessor(), ci.getSema(), false);
        parser->Initialize();
    }

    Clang::~Clang() {
        action->finalizeAction();
    }

    Clang::PartialTranslationUnit& Clang::include(const char* path) {
        using namespace clang;
        Preprocessor& pp = getCI().getPreprocessor();
        assert(pp.isIncrementalProcessingEnabled() && "[Error]: Unexpected clang state");

        std::ostringstream ss{};
        ss << "#include \"" << path << "\"\n";

        const std::string source = ss.str();
        auto memBuf = llvm::WritableMemoryBuffer::getNewUninitMemBuffer(source.length(), path);
        char* memBufStart = memBuf->getBufferStart();
        memcpy(memBufStart, source.data(), source.length());

        SourceManager& SM = getCI().getSourceManager();
        SourceLocation NewLoc = SM.getLocForStartOfFile(SM.getMainFileID());
        FileID FID = SM.createFileID(std::move(memBuf), SrcMgr::C_User, 0, 0, NewLoc);
        if (pp.EnterSourceFile(FID, nullptr, NewLoc))
            throw LLVMException(llvm::make_error<llvm::StringError>("[Error]: Enter source file failed", std::error_code()));

        parse();

        Token AssertTok;
        pp.Lex(AssertTok);
        assert(AssertTok.is(tok::annot_repl_input_end) && "[Error]: Lexer must be EOF when starting incremental parse");
 
        CodeGenerator* codeGen = getCodeGen();
        if (codeGen != nullptr) {
            auto* unitModule = codeGen->ReleaseModule();
            unitModule->setModuleIdentifier(path);
            partialUnitList.front().unitModule.reset(unitModule);
            codeGen->StartModule(DummyFile, unitModule->getContext());
        }
        return partialUnitList.front();
    }

    clang::CodeGenerator* Clang::getCodeGen() noexcept {
        auto* frontendAction = action->getWrapped();
        const bool isCodeGenAction = frontendAction->hasIRSupport();
        if (!isCodeGenAction)
            return nullptr;
        return static_cast<clang::CodeGenAction*>(frontendAction)->getCodeGenerator();
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
            args.push_back(LLVM_RESOURCE_DIR);
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

    void Clang::parse() {
        using namespace clang;
        Sema& sema = ci.getSema();
        llvm::CrashRecoveryContextCleanupRegistrar<Sema> recoverGuard(&sema);
        Sema::GlobalEagerInstantiationScope GlobalInstantiations(sema, true);
        Sema::LocalEagerInstantiationScope LocalInstantiations(sema);

        PartialTranslationUnit unit{};
        ASTContext& astContext = sema.getASTContext();
        astContext.addTranslationUnitDecl();
        unit.unitDecl = astContext.getTranslationUnitDecl();

        // Skip previous eof due to last incremental input.
        const bool skipEOF = parser->getCurToken().is(tok::annot_repl_input_end);
        if (skipEOF) {
            parser->ConsumeAnyToken();
            parser->ExitScope();
            sema.CurContext = nullptr;
            parser->EnterScope(Scope::DeclScope);
            sema.ActOnTranslationUnitScope(parser->getCurScope());
        }

        Parser::DeclGroupPtrTy ADecl;
        Sema::ModuleImportState ImportState;
        bool atEOF = parser->ParseFirstTopLevelDecl(ADecl, ImportState);
        do {
            atEOF = parser->ParseTopLevelDecl(ADecl, ImportState);
        } while(!atEOF);

        DiagnosticsEngine& Diags = getCI().getDiagnostics();
        if (Diags.hasErrorOccurred()) {
            cleanUnit(unit);
            Diags.Reset(true);
            Diags.getClient()->clear();
            throw LLVMException(llvm::make_error<llvm::StringError>("[Error]: Parsing failed", std::error_code()));
        }
        LocalInstantiations.perform();
        GlobalInstantiations.perform();
        partialUnitList.push_front(std::move(unit));
    }

    void Clang::cleanUnit(PartialTranslationUnit& unit) {
        using namespace clang;
        TranslationUnitDecl* unitDecl = unit.unitDecl;
        TranslationUnitDecl* FirstTU = unitDecl->getFirstDecl();
        if (StoredDeclsMap* map = FirstTU->getPrimaryContext()->getLookupPtr()) {
            for (auto ite = map->begin(); ite != map->end(); ++ite) {
                StoredDeclsList& list = ite->second;
                DeclContextLookupResult R = list.getLookupResult();
                for (NamedDecl* D : R)
                    if (D->getTranslationUnitDecl() == unitDecl)
                        list.remove(D);

                if (list.isNull())
                    map->erase(ite);
            }
        }
    }
}
