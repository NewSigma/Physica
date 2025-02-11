/*
 * Copyright 2024-2025 Weibo He.
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
#include "Physica/Core/Exception/LLVMException.h"
#include "Physica/Core/Utils/Container/Array.h"
#include "Physica/Core/Utils/Unreachable.h"
#include "Physica/Python/CXXPtr.h"
#include "Physica/Python/CXXObj.h"
#include "Physica/Python/PhysicaPython.h"
#include "Physica/Python/FFI/FuncInfo.h"

namespace Physica {
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
        const auto pFunc = lookupFunc(clang::GlobalDecl(pDecl->getDestructor(), clang::CXXDtorType::Dtor_Base));
        auto* pDtor = reinterpret_cast<DtorType>(pFunc);
        pDtor(pObj);
        std::free(pObj);
        pObj = nullptr;
    }

    void CXXObj::construct(py::args) {
        assert(pDecl != nullptr);
        assert(pObj == nullptr && "[Error]: Double construct is not allowed");
        using DefaultCtorType = void (*)(void*);
        for (auto ctor : pDecl->ctors()) {
            const bool isDefaultCtor = ctor->getNumCtorInitializers() == 0;
            if (isDefaultCtor) {
                auto* pFunc = lookupFunc(clang::GlobalDecl(ctor, clang::CXXCtorType::Ctor_Base));
                auto* pCtor = reinterpret_cast<DefaultCtorType>(pFunc);
                pObj = allocateObj(pDecl);
                pCtor(pObj);
                break;
            }
        }

        if (pObj == nullptr)
            throw LLVMException("No available default contructor");
    }

    py::object CXXObj::call(const char* rtnTyName, const char* name, py::args args) {
        clang::GlobalDecl funcDecl;
        for (auto pFunc : pDecl->methods()) {
            using namespace clang;
            if (pFunc->getDeclName().getNameKind() != DeclarationName::NameKind::Identifier)
                continue;

            const bool found = pFunc->getName().equals(name);
            if (found) {
                funcDecl = clang::GlobalDecl(pFunc);
                break;
            }
        }
        const auto fn = lookupFunc(std::move(funcDecl));

        auto& pp = PhysicaPython::getInstance();
        const auto rtnType = pp.toCXXType(rtnTyName);
        auto pRtn = rtnType.allocate();

        const size_t numArgs = args.size();
        Array<const ffi_type*> argTypes(numArgs + 1);
        Array<void*> pArgs(numArgs + 1);
        argTypes[0] = &ffi_type_pointer;
        pArgs[0] = &pObj;
        for (size_t i = 1; i <= numArgs; ++i){
            const auto argType = pp.toCXXType(args[i]);
            argTypes[i] = argType.toFFI();
            pArgs[i] = argType.allocate().release();
        }
        FuncInfo info(numArgs + 1, rtnType.toFFI(), argTypes.data());
        ffi_call(const_cast<ffi_cif*>(info.cif()), fn, pRtn.get(), reinterpret_cast<void**>(pArgs.data()));
        
        auto result = rtnType.toPython(pRtn.get());
        for (size_t i = 1; i <= numArgs; ++i)
            std::free(pArgs[i]);
        return result;
    }

    void* CXXObj::allocateObj(const CXXRecordDecl* pDecl) {
        assert(pDecl != nullptr && "[Error]: Invalid param");
        auto& pp = PhysicaPython::getInstance();
        const auto& ctx = pp.getClang().getASTContext();
        const auto type = ctx.getRecordType(pDecl);
        return std::aligned_alloc(ctx.getTypeAlign(type), ctx.getTypeSize(type));
    }

    const CXXObj::CXXRecordDecl* CXXObj::findSpecialization(const ClassTemplateDecl& templateDecl, py::args tparams) {
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
            unreachable();
        }
    }

    CXXObj::ForeignFunc CXXObj::lookupFunc(clang::GlobalDecl decl) {
        auto& pp = PhysicaPython::getInstance();
        auto& codeGen = pp.getClang().getCodeGen();
        const auto symbol = codeGen.GetMangledName(decl);
        auto execAddr = pp.getJIT().lookup(symbol);
        const bool isFound = !bool(execAddr.takeError());
        if (isFound)
            return execAddr.get().toPtr<ForeignFunc>();

        auto& cgm = codeGen.CGM();
        cgm.EmitGlobalDefinition(decl);
        pp.compile(symbol.str().c_str());
        execAddr = pp.getJIT().lookup(symbol);
        return check(std::move(execAddr)).toPtr<ForeignFunc>();
    }
}
