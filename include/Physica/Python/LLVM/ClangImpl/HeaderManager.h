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

#include <unordered_map>
#include "clang/AST/Decl.h"
#include "clang/Lex/PPCallbacks.h"

namespace Physica::Python {
    class Clang;

    class HeaderManager : public clang::PPCallbacks {
        using This = HeaderManager;
        using Base = clang::PPCallbacks;
        using FileID = clang::FileID;
        using StringRef = llvm::StringRef;
        using CharacteristicKind = clang::SrcMgr::CharacteristicKind;
        using HeaderUnitMap = std::unordered_map<std::string, clang::TranslationUnitDecl*>;
        using typename Base::LexedFileChangeReason;

        Clang& compiler;
        HeaderUnitMap loadedHeaders;
    public:
        HeaderManager(Clang& compiler_) : compiler(compiler_) {}
        HeaderManager(const This&) = delete;
        HeaderManager(This&&) noexcept = delete;
        ~HeaderManager() override = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void LexedFileChanged(
                FileID FID,
                LexedFileChangeReason Reason,
                CharacteristicKind FileType,
                FileID PrevFID,
                clang::SourceLocation Loc) override;
        /* Getters */
        [[nodiscard]] const HeaderUnitMap& getLoadedHeaders() const noexcept { return loadedHeaders; }
    private:
        [[nodiscard]] const clang::Preprocessor& getPreprocessor() const noexcept;
    };
}
