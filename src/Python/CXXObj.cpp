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
#include "Physica/Python/FFI/FuncInfo.h"
#include "Physica/Python/CXXPtr.h"
#include "Physica/Python/CXXObj.h"
#include "Physica/Python/PhysicaPython.h"

using namespace Physica;

CXXObj::CXXObj(CXXPtr p, nanobind::args tparams) : pObj(nullptr) {
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

CXXObj::CXXObj(CXXRecordDecl* pDecl_, void* pObj_) noexcept
        : pDecl(pDecl_), pObj(pObj_) {}

CXXObj::CXXObj(CXXObj&& obj) noexcept : pDecl(obj.pDecl), pObj(std::exchange(obj.pObj, nullptr)) {}

CXXObj::~CXXObj() {
    if (pObj == nullptr)
        return;

    using DtorType = void (*)(void*);
    const auto pFunc = lookupFunc(clang::GlobalDecl(pDecl->getDestructor(), clang::CXXDtorType::Dtor_Base));
    auto* pDtor = reinterpret_cast<DtorType>(pFunc);
    pDtor(pObj.get());
    pObj = nullptr;
}

void CXXObj::construct(const nanobind::args&) {
    assert(pDecl != nullptr);
    assert(pObj == nullptr && "[Error]: Double construct is not allowed");
    using DefaultCtorType = void (*)(void*);
    for (auto* ctor : pDecl->ctors()) {
        const bool isDefaultCtor = ctor->getNumCtorInitializers() == 0;
        if (isDefaultCtor) {
            auto* pFunc = lookupFunc(clang::GlobalDecl(ctor, clang::CXXCtorType::Ctor_Base));
            auto* pCtor = reinterpret_cast<DefaultCtorType>(pFunc);
            pObj = allocate(getType());
            pCtor(pObj.get());
            break;
        }
    }

    if (pObj == nullptr)
        throw LLVMException("No available default contructor");
}

nanobind::object CXXObj::call(const char* rtnTyName, const char* name, const nanobind::args& args) {
    clang::GlobalDecl funcDecl;
    for (auto* pFunc : pDecl->methods()) {
        using namespace clang;
        if (pFunc->getDeclName().getNameKind() != DeclarationName::NameKind::Identifier)
            continue;

        const bool found = pFunc->getName().compare(name) == 0;
        if (found) {
            funcDecl = clang::GlobalDecl(pFunc);
            break;
        }
    }
    const auto fn = lookupFunc(funcDecl);

    auto& pp = PhysicaPython::getInstance();
    const auto rtnType = pp.toCXXType(rtnTyName);
    auto pRtn = allocate(rtnType);

    const size_t numArgs = args.size();
    Array<const ffi_type*> argTypes(numArgs + 1);
    Array<Ptr> pArgs(numArgs + 1);
    Array<void*> pRawArgs(numArgs + 1);
    argTypes[0] = &ffi_type_pointer;
    pRawArgs[0] = pObj.get();
    for (size_t i = 1; i <= numArgs; ++i) {
        const auto argType = pp.toCXXType(args[i]);
        argTypes[i] = argType.toFFI();
        pArgs[i] = allocate(argType);
        pRawArgs[i] = pArgs[i].get();
    }
    FuncInfo info(numArgs + 1, rtnType.toFFI(), argTypes.data());
    ffi_call(const_cast<ffi_cif*>(info.cif()), fn, pRtn.get(), pRawArgs.data());
    return rtnType.toPython(pRtn.get());
}

auto CXXObj::allocate(CXXType ty) -> Ptr {
    return Ptr(::operator new(ty.getSize(), std::align_val_t(ty.getAlign())));
}

CXXObj::CXXRecordDecl* CXXObj::findSpecialization(const ClassTemplateDecl& templateDecl, nanobind::args tparams) {
    using namespace clang;
    const size_t numArgs = tparams.size();
    for (ClassTemplateSpecializationDecl* pSpec : templateDecl.specializations()) {
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

bool CXXObj::matchParamT(const nanobind::handle& pyarg, const clang::TemplateArgument& targ) {
    using Kind = clang::TemplateArgument::ArgKind;
    switch (targ.getKind()) {
    case Kind::Null:
    case Kind::Type:
    case Kind::Declaration:
    case Kind::NullPtr:
        return false;
    case Kind::Integral:
        return targ.getAsIntegral() == int64_t(nanobind::int_(pyarg));
    case Kind::Template:
    case Kind::TemplateExpansion:
    case Kind::Expression:
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
