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
#include <iostream>
#include "Physica/Python/LLVM/ClangImpl/HeaderManager.h"
#include "Physica/Python/LLVM/Clang.h"

namespace Physica::Python {
    void HeaderManager::LexedFileChanged(
            FileID,
            LexedFileChangeReason Reason,
            CharacteristicKind FileType,
            FileID PrevFID,
            clang::SourceLocation) {
        const bool isUserHeader = FileType == CharacteristicKind::C_User;
        const bool isDone = Reason == LexedFileChangeReason::ExitFile;
        if (!isUserHeader || !isDone)
            return;

        const auto* pFile = getPreprocessor().getSourceManager().getFileEntryForID(PrevFID);
        if (pFile == nullptr)
            return;

        constexpr const char IncludeDir[] = "include/Physica/";
        const StringRef fileName = pFile->getName();
        const size_t pos = fileName.find(IncludeDir);
        const bool isNotPhysicaHeader = pos == StringRef::npos;
        if (isNotPhysicaHeader)
            return;
        
        const StringRef includeName = fileName.substr(pos + sizeof(IncludeDir) - 1);
        const bool isImpl = includeName.find("Impl") != StringRef::npos;
        if (isImpl)
            return;

        StringRef className;
        /* Make class name */ {
            const size_t end = includeName.rfind('.');
            assert(end != StringRef::npos && "[Error]: Unexpected includeName");
            className = includeName.slice(includeName.rfind('/') + 1, end);
        }
        auto& ctx = compiler.getCI().getASTContext();
        compiler.getCI().getASTConsumer().HandleTranslationUnit(ctx);
        auto* pClass = findClassDecl(*ctx.getTranslationUnitDecl(), className);
        loadedHeaders[includeName.str()] = pClass;
    }

    const clang::Preprocessor& HeaderManager::getPreprocessor() const noexcept {
        return compiler.getCI().getPreprocessor();
    }

    HeaderManager::ClassDecl* HeaderManager::findClassDecl(DeclContext& tree, StringRef className) {
        assert(!className.empty() && "[Error]: Invalid class name");
        using namespace clang;
        for (clang::Decl* pDecl : tree.decls()) {
            if (!llvm::isa<NamedDecl>(*pDecl))
                continue;

            NamedDecl& decl = static_cast<NamedDecl&>(*pDecl);
            if (!decl.getDeclName().isIdentifier())
                continue;

            const auto declName = decl.getName();
            if (llvm::isa<NamespaceDecl>(decl)) {
                const bool isNotPhysicaNamespace = declName.empty() || std::islower(declName.front());
                if (isNotPhysicaNamespace)
                    continue;

                auto* result = findClassDecl(static_cast<NamespaceDecl&>(decl), className);
                if (result != nullptr)
                    return result;
            }
            else if (llvm::isa<CXXRecordDecl>(decl) || llvm::isa<ClassTemplateDecl>(decl)) {
                const bool found = className.compare(declName) == 0;
                if (found)
                    return &decl;
            }
        }
        return nullptr;
    }
}
