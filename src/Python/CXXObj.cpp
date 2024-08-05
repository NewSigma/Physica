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
#include <clang/AST/GlobalDecl.h>
#include <Physica/Python/LLVM/LLVM.h>
#include <Physica/Python/CXXObj.h>
#include <Physica/Python/PhysicaPython.h>
#include <Physica/Python/Exception/LLVMException.h>
#include <Physica/Python/FFI/FuncInfo.h>
#include <Physica/Utils/Container/Array/Array.h>

namespace Physica::Python {
    CXXObj::CXXObj(const CXXRecordDecl* pDecl_, void* pObj_) noexcept : pDecl(pDecl_), pObj(pObj_) {}

    CXXObj::CXXObj(CXXObj&& obj) noexcept : pDecl(obj.pDecl), pObj(obj.pObj) {
        obj.pObj = nullptr;
    }

    CXXObj::~CXXObj() {
        if (pObj == nullptr)
            return;

        using DtorType = void (*)(void*);
        auto& pp = PhysicaPython::getInstance();
        Clang& clang = pp.getClang();
        auto dtor = clang::GlobalDecl(pDecl->getDestructor(), clang::CXXDtorType::Dtor_Base);
        const std::string symbol = clang.getCodeGen().GetMangledName(std::move(dtor)).str();
        const auto pDtor = llvmCheck(pp.getJIT().lookup(symbol));
        pDtor.toPtr<DtorType>()(pObj);
        std::free(pObj);
        pObj = nullptr;
    }

    py::object CXXObj::call(const char* rtnTyName, const char* name, py::args args) {
        clang::GlobalDecl funcDecl;
        for (auto pFunc : pDecl->methods()) {
            const bool found = pFunc->getName() == llvm::StringRef(name);
            if (found)
                funcDecl = clang::GlobalDecl(pFunc);
        }
        using ForeignFunc = void (*)();
        auto& pp = PhysicaPython::getInstance();
        Clang& clang = pp.getClang();
        const std::string symbol = clang.getCodeGen().GetMangledName(std::move(funcDecl)).str();
        const auto fn = llvmCheck(pp.getJIT().lookup(symbol)).toPtr<ForeignFunc>();

        const auto rtnType = pp.toCXXType(rtnTyName);
        auto pRtn = rtnType.allocate();

        const size_t numArgs = args.size();
        Utils::Array<const ffi_type*> argTypes(numArgs + 1);
        Utils::Array<Core::PlainPtr> pArgs(numArgs + 1);
        argTypes[0] = &ffi_type_pointer;
        pArgs[0] = Core::PlainPtr(&pObj);
        for (size_t i = 1; i <= numArgs; ++i){
            const auto argType = pp.toCXXType(args[i]);
            argTypes[i] = argType.toFFI();
            pArgs[i] = argType.allocate();
        }
        FuncInfo info(numArgs + 1, rtnType.toFFI(), argTypes.data());
        ffi_call(const_cast<ffi_cif*>(info.cif()), fn, pRtn.get(), reinterpret_cast<void**>(pArgs.data()));
        pArgs[0] = nullptr;
        return rtnType.toPython(std::move(pRtn));
    }

    CXXObj CXXObj::create(const CXXRecordDecl* pDecl) {
        assert(pDecl != nullptr); 
        using DefaultCtorType = void (*)(void*);
        auto& pp = PhysicaPython::getInstance();
        Clang& clang = pp.getClang();
        for (auto ctor : pDecl->ctors()) {
            const bool isDefaultCtor = ctor->getNumCtorInitializers() == 0;
            if (isDefaultCtor) {
                auto defaultCtor = clang::GlobalDecl(ctor, clang::CXXCtorType::Ctor_Base);
                const std::string symbol = clang.getCodeGen().GetMangledName(std::move(defaultCtor)).str();
                const auto pCtor = llvmCheck(pp.getJIT().lookup(symbol));
                void* pObj = allocateObj(pDecl);
                pCtor.toPtr<DefaultCtorType>()(pObj);
                return CXXObj(pDecl, pObj);
            }
        }
        throw LLVMException("No available default contructor");
    }

    void* CXXObj::allocateObj(const CXXRecordDecl* pDecl) {
        assert(pDecl != nullptr && "[Error]: Invalid param");
        const auto& pp = PhysicaPython::getInstance();
        const auto& ctx = pp.getClang().getASTContext();
        const auto type = ctx.getRecordType(pDecl);
        return std::aligned_alloc(ctx.getTypeAlign(type), ctx.getTypeSize(type));
    }
}
