/*
 * Copyright 2024-2026 Weibo He.
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
#include "clang/AST/DeclContextInternals.h"
#include "clang/Basic/Diagnostic.h"
#include "clang/CodeGen/CodeGenAction.h"
#include "clang/Driver/Driver.h"
#include "clang/Driver/Compilation.h"
#include "clang/Driver/Tool.h"
#include "clang/Frontend/TextDiagnosticPrinter.h"
#include "clang/Lex/PreprocessorOptions.h"
#include "clang/Lex/Preprocessor.h"
#include "clang/Sema/Sema.h"
#include "llvm/TargetParser/Host.h"
#include "llvm/IR/Module.h"
#include "llvm/Support/CrashRecoveryContext.h"
#include "Physica/Core/Exception/LLVMException.h"
#include "Physica/Python/LLVM/Clang.h"
#include "Physica/Python/LLVM/ClangImpl/HeaderManager.h"
#include "Physica/Macro.h"

using namespace Physica;

Clang::Clang(std::filesystem::path root_, LLVM& llvm)
        : root(std::move(root_)) {
    makeInvocation();
    initOptions();
    Base::createTarget();
    Base::createDiagnostics();
    Base::LoadRequestedPlugins();

    action = std::make_unique<IncrementalAction>(*this, llvm.getContext());
    Base::ExecuteAction(*action);

    auto& pp = Base::getPreprocessor();
    pHeaderManager = new HeaderManager(*this);
    pp.addPPCallbacks(std::unique_ptr<clang::PPCallbacks>(pHeaderManager));

    parser = std::make_unique<Parser>(pp, Base::getSema(), false);
    parser->Initialize();
}

const clang::NamedDecl* Clang::include(const char* includeName) {
    const auto& headers = getHeaderManager().getLoadedHeaders();
    auto ite = headers.find(std::string(includeName));
    const bool isNotFound = ite == headers.cend();
    if (isNotFound) {
        using namespace clang;
        std::unique_ptr<llvm::WritableMemoryBuffer> memBuf;
        /* Make memBuf */ {
            const std::string source = std::format("#include \"Physica/{}\"\n", includeName);
            memBuf = llvm::WritableMemoryBuffer::getNewUninitMemBuffer(source.length(), includeName);
            memcpy(memBuf->getBufferStart(), source.data(), source.length());
        }
        SourceManager& manager = Base::getSourceManager();
        SourceLocation loc = manager.getLocForStartOfFile(manager.getMainFileID());
        FileID fid = manager.createFileID(std::move(memBuf), SrcMgr::C_User, 0, 0, loc);

        Preprocessor& pp = Base::getPreprocessor();
        assert(pp.isIncrementalProcessingEnabled() && "[Error]: Unexpected clang state");
        const bool failed = pp.EnterSourceFile(fid, nullptr, loc);
        if (failed)
            throw LLVMException("[Error]: Enter source file failed");
        parse();
        ite = headers.find(std::string(includeName));
    }
    assert(ite != headers.cend());
    return (*ite).second;
}

Clang::PartialTranslationUnit& Clang::compile(const char* moduleName) {
    using namespace clang;
    Sema& sema = Base::getSema();
    llvm::CrashRecoveryContextCleanupRegistrar<Sema> recoverGuard(&sema);

    auto& ctx = getASTContext();
    ctx.addTranslationUnitDecl();
    parse();

    PartialTranslationUnit ptu{};
    ptu.unitDecl = ctx.getTranslationUnitDecl();
    /* Make module */ {
        CodeGenerator& codeGen = getCodeGen();
        ptu.unitModule = codeGen.ReleaseModule();
        ptu.unitModule->setModuleIdentifier(moduleName);
        codeGen.StartModule(DummyFile, ptu.unitModule->getContext());
    }
    partialUnitList.push_front(std::move(ptu));
    return partialUnitList.front();
}

clang::CodeGenerator& Clang::getCodeGen() noexcept {
    auto* frontendAction = action->getWrapped();
    assert(frontendAction->hasIRSupport() && "[Error]: No CodeGen because it is not CodeGenAction");
    return *static_cast<clang::CodeGenAction*>(frontendAction)->getCodeGenerator();
}

void Clang::makeInvocation() {
    const std::string exec = llvm::sys::fs::getMainExecutable(nullptr, nullptr);
    std::filesystem::path root_inc;
    std::vector<const char*> args{};
    /* Make args */ {
        args.push_back(exec.c_str());
        constexpr bool NoImpl = false;
        if constexpr (NoImpl) {
            args.push_back("-xcuda");
            args.push_back("--cuda-host-only");
            args.push_back(PHYSICA_CUDA_ARCHITECTURES);
        }
        else
            args.push_back("-xc++");
        args.push_back("-std=c++23");
        args.push_back("-Xclang");
        args.push_back("-fincremental-extensions");
        args.push_back("--compile");
        args.push_back(DummyFile);
        args.push_back("-resource-dir");
        args.push_back(LLVM_RESOURCE_DIR);
        {
            root_inc = root;
            root_inc.append("include");
            args.push_back("-I");
            args.push_back(root_inc.c_str());
        }
        if constexpr (HasHDF5()) {
            args.push_back("-I");
            args.push_back(HDF5_INCLUDE_DIR);
        }
    }

    using namespace clang;
    {
        auto pDiagOpts = CreateAndPopulateDiagOpts(args);
        auto* pDiagBuffer = new TextDiagnosticPrinter(llvm::outs(), *pDiagOpts);
        auto* pDiag = new DiagnosticsEngine(new DiagnosticIDs(), *pDiagOpts.release(), pDiagBuffer);
        Base::setDiagnostics(pDiag);
    }
    auto llvmDirver = driver::Driver("", llvm::sys::getProcessTriple(), Base::getDiagnostics());
    std::unique_ptr<driver::Compilation> pCompile;
    /* Make compilation */ {
        llvmDirver.setCheckInputsExist(false);
        pCompile.reset(llvmDirver.BuildCompilation(args));
        if (pCompile->getArgs().hasArg(options::OPT_v))
            pCompile->getJobs().Print(llvm::errs(), "\n", false);
    }

    const auto& jobs = pCompile->getJobs();
    if ((jobs.size() == 0) || !isa<driver::Command>(*jobs.begin()))
        throw LLVMException("[Error]: Unable to create a driver job");

    const driver::Command* cmd = cast<driver::Command>(&(*jobs.begin()));
    if (strcmp("clang", cmd->getCreator().getName()) != 0)
        throw LLVMException("[Error]: Driver initialization failed");

    const bool success = CompilerInvocation::CreateFromArgs(Base::getInvocation(), cmd->getArguments(), Base::getDiagnostics());
    if (!success)
        throw LLVMException("[Error]: Failed to create compiler invocation");
}

void Clang::initOptions() {
    using namespace clang;
    /* Preprocessor */ {
        auto* memBuffer = llvm::MemoryBuffer::getMemBuffer("").release();
        Base::getPreprocessorOpts().addRemappedFile(DummyFile, memBuffer);
    }
    /* CodeGen */ {
        auto& options = Base::getCodeGenOpts();
        options.ClearASTBeforeBackend = false;
        options.DisableFree = false;
    }
    Base::getFrontendOpts().DisableFree = false;
}

void Clang::parse() {
    using namespace clang;
    Sema& sema = Base::getSema();
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
        auto& consumer = Base::getASTConsumer();
        Parser::DeclGroupPtrTy pDeclGroup;
        Sema::ModuleImportState state{};
        bool isDone = parser->ParseFirstTopLevelDecl(pDeclGroup, state);
        while (!isDone) {
            if (pDeclGroup) {
                auto declGroup = pDeclGroup.get();
                if (!consumer.HandleTopLevelDecl(declGroup))
                    throw LLVMException("[Error]: Consumer rejected the decl");
            }
            isDone = parser->ParseTopLevelDecl(pDeclGroup, state);
        }
    }

    auto& diag = Base::getDiagnostics();
    if (diag.hasErrorOccurred()) {
        cleanLastUnit();
        diag.Reset(true);
        diag.getClient()->clear();
        throw LLVMException("[Error]: Parsing failed");
    }
    Token check{};
    Base::getPreprocessor().Lex(check);
    assert(check.is(tok::annot_repl_input_end) && "[Error]: Lexer must be EOF when starting incremental parse");
}

void Clang::cleanLastUnit() noexcept {
    using namespace clang;
    auto* unitDecl = Base::getSema().getASTContext().getTranslationUnitDecl();
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
