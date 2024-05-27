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

        auto& ctx = compiler.getCI().getASTContext();
        compiler.getCI().getASTConsumer().HandleTranslationUnit(ctx);
        loadedHeaders[includeName.str()] = ctx.getTranslationUnitDecl();
        ctx.addTranslationUnitDecl();
    }

    const clang::Preprocessor& HeaderManager::getPreprocessor() const noexcept {
        return compiler.getCI().getPreprocessor();
    }
}
