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
#include "Physica/Python/LLVM/LLVM.h"
#include "Physica/Python/LLVM/Clang.h"
#include "Physica/Python/LLVM/ClangImpl/HeaderManager.h"
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
        makeInvocation();
        makeOptions();
        ci.createDiagnostics();
        ci.LoadRequestedPlugins();

        action = std::make_unique<IncrementalAction>(ci, LLVM::getInstance().getContext());
        ci.ExecuteAction(*action);

        auto& pp = ci.getPreprocessor();
        pp.addPPCallbacks(std::unique_ptr<PPCallbacks>(new HeaderManager(*this)));
        parser = std::make_unique<Parser>(pp, ci.getSema(), false);
        parser->Initialize();
    }

    Clang::~Clang() {
        action->finalizeAction();
    }

    const typename Clang::TranslationUnitDecl* Clang::include(const char* includeName) {
        const auto& headers = getHeaderManager().getLoadedHeaders();
        const auto ite = headers.find(includeName);
        if (ite != headers.cend())
            return ite->second;

        using namespace clang;
        std::unique_ptr<llvm::WritableMemoryBuffer> memBuf;
        /* Make memBuf */ {
            std::ostringstream ss{};
            ss << "#include \"Physica/" << includeName << "\"\n";
            const std::string source = ss.str();
            memBuf = llvm::WritableMemoryBuffer::getNewUninitMemBuffer(source.length(), includeName);
            memcpy(memBuf->getBufferStart(), source.data(), source.length());
        }
        SourceManager& manager = getCI().getSourceManager();
        SourceLocation loc = manager.getLocForStartOfFile(manager.getMainFileID());
        FileID fid = manager.createFileID(std::move(memBuf), SrcMgr::C_User, 0, 0, loc);

        Preprocessor& pp = getCI().getPreprocessor();
        assert(pp.isIncrementalProcessingEnabled() && "[Error]: Unexpected clang state");
        const bool failed = pp.EnterSourceFile(std::move(fid), nullptr, std::move(loc));
        if (failed)
            throw LLVMException("[Error]: Enter source file failed");
        return makePTU(includeName).unitDecl;
    }

    Clang::PartialTranslationUnit& Clang::makePTU(const char* moduleName) {
        using namespace clang;
        Sema& sema = ci.getSema();
        llvm::CrashRecoveryContextCleanupRegistrar<Sema> recoverGuard(&sema);
        Sema::GlobalEagerInstantiationScope GlobalInstantiations(sema, true);
        Sema::LocalEagerInstantiationScope LocalInstantiations(sema);

        auto& ctx = sema.getASTContext();
        ctx.addTranslationUnitDecl();
        parse();

        LocalInstantiations.perform();
        GlobalInstantiations.perform();
        ci.getASTConsumer().HandleTranslationUnit(ctx);

        PartialTranslationUnit unit{};
        unit.unitDecl = ctx.getTranslationUnitDecl();
        CodeGenerator* codeGen = getCodeGen();
        if (codeGen != nullptr) {
            auto* unitModule = codeGen->ReleaseModule();
            unitModule->setModuleIdentifier(moduleName);
            unit.unitModule.reset(unitModule);
            codeGen->StartModule(DummyFile, unitModule->getContext());
        }
        partialUnitList.push_front(std::move(unit));
        return partialUnitList.front();
    }

    clang::CodeGenerator* Clang::getCodeGen() noexcept {
        auto* frontendAction = action->getWrapped();
        const bool isCodeGenAction = frontendAction->hasIRSupport();
        if (!isCodeGenAction)
            return nullptr;
        return static_cast<clang::CodeGenAction*>(frontendAction)->getCodeGenerator();
    }

    const HeaderManager& Clang::getHeaderManager() const noexcept {
        return static_cast<HeaderManager&>(*ci.getPreprocessor().getPPCallbacks());
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
            args.push_back("-resource-dir");
            args.push_back(LLVM_RESOURCE_DIR);
            args.push_back("-I");
            args.push_back(PHYSICA_INCLUDE_DIR);
            if constexpr (IsHDF5Enabled()) {
                args.push_back("-I");
                args.push_back(HDF5_INCLUDE_DIR);
            }
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

        auto* diagOptions = new DiagnosticOptions();
        auto* diagBuffer1 = new TextDiagnosticPrinter(llvm::outs(), diagOptions);
        auto* diag = new DiagnosticsEngine(new DiagnosticIDs(), diagOptions, diagBuffer1, true);
        const bool success = CompilerInvocation::CreateFromArgs(ci.getInvocation(), cmd->getArguments(), *diag);
        if (!success)
            throw LLVMException(llvm::createStringError(llvm::errc::not_supported, "[Error]: Failed to create compiler invocation"));
    }

    void Clang::makeOptions() {
        using namespace clang;
        auto* memBuffer = llvm::MemoryBuffer::getMemBuffer("").release();
        ci.getPreprocessorOpts().addRemappedFile(DummyFile, memBuffer);
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
        // Skip previous eof due to last incremental input.
        const bool skipEOF = parser->getCurToken().is(tok::annot_repl_input_end);
        if (skipEOF) {
            parser->ConsumeAnyToken();
            parser->ExitScope();
            sema.CurContext = nullptr;
            parser->EnterScope(Scope::DeclScope);
            sema.ActOnTranslationUnitScope(parser->getCurScope());
        }

        {
            auto& consumer = ci.getASTConsumer();
            Parser::DeclGroupPtrTy declGroup;
            Sema::ModuleImportState state;
            bool isDone = parser->ParseFirstTopLevelDecl(declGroup, state);
            while (!isDone) {
                if (declGroup && !consumer.HandleTopLevelDecl(declGroup.get()))
                    throw LLVMException("[Error]: Consumer rejected the decl");
                isDone = parser->ParseTopLevelDecl(declGroup, state);
            }
        }
        DiagnosticsEngine& Diags = getCI().getDiagnostics();
        if (Diags.hasErrorOccurred()) {
            cleanLastUnit();
            Diags.Reset(true);
            Diags.getClient()->clear();
            throw LLVMException("[Error]: Parsing failed");
        }
        Token check;
        getCI().getPreprocessor().Lex(check);
        assert(check.is(tok::annot_repl_input_end) && "[Error]: Lexer must be EOF when starting incremental parse");
    }

    void Clang::cleanLastUnit() noexcept {
        using namespace clang;
        auto* unitDecl = ci.getSema().getASTContext().getTranslationUnitDecl();
        TranslationUnitDecl* FirstTU = unitDecl->getFirstDecl();
        if (StoredDeclsMap* map = FirstTU->getPrimaryContext()->getLookupPtr()) {
            for (auto ite = map->begin(); ite != map->end(); ++ite) {
                StoredDeclsList& list = ite->second;
                DeclContextLookupResult R = list.getLookupResult();
                for (NamedDecl* D : R)
                    if (D != nullptr && D->getTranslationUnitDecl() == unitDecl)
                        list.remove(D);

                if (list.isNull())
                    map->erase(ite);
            }
        }
    }
}
