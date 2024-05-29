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
#include "clang/AST/GlobalDecl.h"
#include "Physica/Python/LLVM/LLVM.h"
#include "Physica/Python/CXXObj.h"
#include "Physica/Python/PhysicaPython.h"
#include "Physica/Python/Exception/LLVMException.h"

namespace Physica::Python {
    CXXObj::CXXObj(const CXXRecordDecl* pDecl_, void* pObj_) : pDecl(pDecl_), pObj(pObj_) {}

    CXXObj::CXXObj(CXXObj&& obj) noexcept : pDecl(obj.pDecl), pObj(obj.pObj) {
        obj.pObj = nullptr;
    }

    CXXObj::~CXXObj() {
        if (pObj == nullptr)
            return;

        using DestructorType = void (*)(void*);
        auto& pp = PhysicaPython::getInstance();
        Clang& clang = pp.getClang();
        auto destructor = clang::GlobalDecl(pDecl->getDestructor(), clang::CXXDtorType::Dtor_Base);
        const std::string symbol = clang.getCodeGen().GetMangledName(std::move(destructor)).str();
        const auto pDel = llvmCheck(pp.getJIT().lookup(symbol));
        pDel.toPtr<DestructorType>()(pObj);
        pObj = nullptr;
    }

    py::object CXXObj::call(py::handle type, const char* name) {
        const std::string type_str = py::str(type);
        printf("%s|%s", type_str.c_str(), name);
        return py::none();
    }
}
