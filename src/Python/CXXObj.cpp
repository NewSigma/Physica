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
#include <clang/AST/GlobalDecl.h>
#include <Physica/Python/LLVM/LLVM.h>
#include <Physica/Python/CXXPtr.h>
#include <Physica/Python/CXXObj.h>
#include <Physica/Python/PhysicaPython.h>
#include <Physica/Python/Exception/LLVMException.h>
#include <Physica/Python/FFI/FuncInfo.h>
#include <Physica/Utils/Container/Array/Array.h>
#include <Physica/Utils/Unreachable.h>

namespace Physica::Python {
    CXXObj::CXXObj(CXXPtr p, py::args tparams) : pObj(nullptr) {
        using namespace clang;
        const bool isPlainClass = tparams.size() == 0;
        if (isPlainClass) {
            assert(llvm::isa<CXXRecordDecl>(*p.cast<Decl>()) && "[Error]: This is not a plain class");
            pDecl = p.cast<CXXRecordDecl>();
        }
        else {
            assert(llvm::isa<ClassTemplateDecl>(*p.cast<Decl>()) && "[Error]: This is not a template class");
            pDecl = findSpecialization(*p.cast<ClassTemplateDecl>(), std::move(tparams));
            if (pDecl == nullptr)
                throw LLVMException("[Error]: Cannot find specialization");
        }
    }

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

    void CXXObj::construct(py::args) {
        assert(pDecl != nullptr);
        assert(pObj == nullptr && "[Error]: Double construct is not allowed");
        using DefaultCtorType = void (*)(void*);
        auto& pp = PhysicaPython::getInstance();
        Clang& clang = pp.getClang();
        for (auto ctor : pDecl->ctors()) {
            const bool isDefaultCtor = ctor->getNumCtorInitializers() == 0;
            if (isDefaultCtor) {
                auto defaultCtor = clang::GlobalDecl(ctor, clang::CXXCtorType::Ctor_Base);
                const std::string symbol = clang.getCodeGen().GetMangledName(std::move(defaultCtor)).str();
                const auto pCtor = llvmCheck(pp.getJIT().lookup(symbol));
                pObj = allocateObj(pDecl);
                pCtor.toPtr<DefaultCtorType>()(pObj);
                break;
            }
        }

        if (pObj == nullptr)
            throw LLVMException("No available default contructor");
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

    void* CXXObj::allocateObj(const CXXRecordDecl* pDecl) {
        assert(pDecl != nullptr && "[Error]: Invalid param");
        const auto& pp = PhysicaPython::getInstance();
        const auto& ctx = pp.getClang().getASTContext();
        const auto type = ctx.getRecordType(pDecl);
        return std::aligned_alloc(ctx.getTypeAlign(type), ctx.getTypeSize(type));
    }

    const typename CXXObj::CXXRecordDecl* CXXObj::findSpecialization(const ClassTemplateDecl& templateDecl, py::args tparams) {
        using namespace clang;
        const size_t numArgs = tparams.size();
        for (const ClassTemplateSpecializationDecl* pSpec : templateDecl.specializations()) {
            const auto& specArgs = pSpec->getTemplateArgs();
            if (numArgs != specArgs.size())
                continue;
            
            bool match = true;
            for (size_t i = 0; i < numArgs && match; ++i)
                match &= matchParamT(tparams[i], specArgs[i]);

            if (match)
                return pSpec;
        }
        return nullptr;
    }

    bool CXXObj::matchParamT(const py::handle& pyarg, const clang::TemplateArgument& targ) {
        using Kind = clang::TemplateArgument::ArgKind;
        switch (targ.getKind()) {
        case Kind::Null:
            return false;
        case Kind::Type:
            return false;
        case Kind::Declaration:
            return false;
        case Kind::NullPtr:
            return false;
        case Kind::Integral:
            return targ.getAsIntegral() == int64_t(pyarg.cast<py::int_>());
        case Kind::Template:
            return false;
        case Kind::TemplateExpansion:
            return false;
        case Kind::Expression:
            return false;
        case Kind::Pack:
            return false;
        default:
            assert(false && "[Error]: Invalid tparam type");
            Utils::unreachable();
        }
    }
}
