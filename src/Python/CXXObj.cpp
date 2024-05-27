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
    CXXObj::CXXObj(CXXObj&& obj) noexcept : decl(obj.decl), pObj(obj.pObj) {
        obj.pObj = nullptr;
    }

    CXXObj::~CXXObj() {
        //using DestructorType = void (*)(void*);
        if (pObj == nullptr)
            return;

        //pDel(pObj);
        pObj = nullptr;
    }

    py::object CXXObj::call(py::handle type, const char* name) {
        const std::string type_str = py::str(type);
        printf("%s|%s", type_str.c_str(), name);
        return py::none();
    }

    clang::FunctionDecl* CXXObj::makeDestructorAST() {
        using namespace clang;
        auto& ci = PhysicaPython::getInstance().getClang().getCI();
        auto& pp = ci.getPreprocessor();
        auto& idTable = pp.getIdentifierTable();
        ASTContext& ctx = ci.getSema().getASTContext();
        ctx.addTranslationUnitDecl();
        auto* unitDecl = ctx.getTranslationUnitDecl();

        const QualType paramType = ctx.VoidPtrTy;
        const QualType funcType = ctx.getFunctionType(ctx.VoidPtrTy, {paramType}, FunctionProtoType::ExtProtoInfo());
        auto* result = FunctionDecl::Create(ctx, unitDecl, {}, {}, &idTable.getOwn("Del")
                                        , funcType, ctx.CreateTypeSourceInfo(funcType), clang::SC_None);
        unitDecl->addDecl(result);

        auto* pObjParm = ParmVarDecl::Create(ctx, unitDecl, {}, {}, &idTable.getOwn("pObj"),
                                paramType, ctx.CreateTypeSourceInfo(paramType),
                                clang::SC_None, nullptr);
        result->setParams({pObjParm});
/*
       `-CXXReinterpretCastExpr 0x5654a4aff9a8 <col:5, col:30> 'T *' reinterpret_cast<T *> <BitCast>
          `-ImplicitCastExpr 0x5654a4aff990 <col:26> 'void *' <LValueToRValue> part_of_explicit_cast
            `-DeclRefExpr 0x5654a4aff930 <col:26> 'void *' lvalue ParmVar 0x5654a4aff778 'pObj' 'void *'
*/
/*
        auto* body = CompoundStmt::CreateEmpty(ctx, 1, false);
        auto* refParm = DeclRefExpr::Create(ctx, {}, {}, pObjParm, false, SourceLocation{}, paramType, ExprValueKind::VK_LValue);
        auto* valueCastExpr = ImplicitCastExpr::Create(ctx, ctx.VoidPtrTy, CastKind::CK_LValueToRValue, refParm, nullptr, ExprValueKind::VK_PRValue, {});
        auto* ptrCastExpr = CXXReinterpretCastExpr::Create(ctx, QuanType(), ExprValueKind::VK_PRValue, CastKind::CK_BitCast,
                valueCastExpr, nullptr, nullptr, {}, {}, {});
        auto* member = MemberExpr::Create();
        auto* caller = CXXMemberCallExpr::Create(ctx, member, {}, ctx.VoidTy, ExprValueKind::VK_PRValue, {}, {});
        (*body->body_begin()) = caller;

        result->setBody(body);
*/
        if (!ci.getASTConsumer().HandleTopLevelDecl(DeclGroupRef(result))) [[unlikely]]
            throw Physica::Python::LLVMException("[Error]: Consumer rejected the decl");
        return result;
    }
}
