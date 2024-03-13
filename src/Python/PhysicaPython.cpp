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
#include <pybind11/pybind11.h>
#include "llvm/IR/IRBuilder.h"
#include "clang/AST/GlobalDecl.h"
#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Python/PhysicaPython.h"
#include "Physica/Python/LLVM/LLVM.h"
#include "Physica/Python/LLVM/Clang.h"
#include "Physica/Python/LLVM/Executor.h"
#include "Physica/Python/Exception/LLVMException.h"

using namespace Physica::Python;
namespace py = pybind11;

clang::FunctionDecl* makeAST(clang::CompilerInstance& ci);

PYBIND11_MODULE(PhysicaPython, m) {
    m.doc() = "PhysicaPython is a python interface to Physica";
    py::register_exception<LLVMException>(m, "LLVMException");

    using ScalarOption = Physica::Core::ScalarOption;
    py::enum_<ScalarOption>(m, "ScalarOption", py::arithmetic())
        .value("Float", ScalarOption::Float)
        .value("Double", ScalarOption::Double)
        .value("MultiPrecision", ScalarOption::MultiPrecision)
        .export_values();

    m.def("runKernel", []() {
        Clang temp{};
        std::string funcSym;
        auto& plu = temp.makePLU("test", [&]() {
            auto* funcDecl = makeAST(temp.getCI());
            funcSym = temp.getCodeGen()->GetMangledName(clang::GlobalDecl(funcDecl)).str();
        });

        Executor exec(temp.getCI().getTarget());
        auto& jit = exec.getJIT();
        llvmCheck(jit.addIRModule(llvm::orc::ThreadSafeModule(std::move(plu.unitModule), LLVM::getInstance().getThreadSafeContext())));
        auto pAdd1 = llvmCheck(jit.lookup(funcSym));
        int (*add1)(int) = pAdd1.toPtr<int(int)>();
        int result = add1(42);
        std::cout << "add1(42) = " << result << "\n";
    });
}

clang::FunctionDecl* makeAST(clang::CompilerInstance& ci) {
    using namespace clang;
    auto& pp = ci.getPreprocessor();
    auto& idTable = pp.getIdentifierTable();
    auto& astContext = ci.getSema().getASTContext();
    astContext.addTranslationUnitDecl();
    auto* unitDecl = astContext.getTranslationUnitDecl();

    const QualType paramType = astContext.IntTy;
    const QualType funcType = astContext.getFunctionType(paramType, {paramType}, FunctionProtoType::ExtProtoInfo());
    auto* add1 = FunctionDecl::Create(astContext, unitDecl, {}, {}, &idTable.getOwn("add1")
                                    , funcType, astContext.CreateTypeSourceInfo(funcType), clang::SC_None);
    unitDecl->addDecl(add1);

    auto* i = ParmVarDecl::Create(astContext, unitDecl, {}, {}, &idTable.getOwn("i"),
                             paramType, astContext.CreateTypeSourceInfo(paramType),
                             clang::SC_None, nullptr);
    add1->setParams({i});

    auto* body = CompoundStmt::CreateEmpty(astContext, 1, false);
    add1->setBody(body);

    auto* rtnStmt = ReturnStmt::Create(astContext, {}, nullptr, nullptr);
    {
        auto* binaryOp = BinaryOperator::CreateEmpty(astContext, false);
        {
            auto* op = DeclRefExpr::Create(astContext, {}, {}, i, false, SourceLocation{}, paramType, ExprValueKind::VK_LValue);
            auto* lhs = ImplicitCastExpr::Create(astContext, astContext.IntTy, CastKind::CK_LValueToRValue, op, nullptr, ExprValueKind::VK_PRValue, {});
            binaryOp->setLHS(lhs);
        }
        auto* rhs = IntegerLiteral::Create(astContext, llvm::APInt(sizeof(int) * 8, 1, true), astContext.IntTy, {});
        binaryOp->setRHS(rhs);
        binaryOp->setOpcode(BinaryOperator::Opcode::BO_Add);
        binaryOp->setType(astContext.IntTy);
        binaryOp->setValueKind(ExprValueKind::VK_PRValue);
        binaryOp->setObjectKind(ExprObjectKind::OK_Ordinary);
        rtnStmt->setRetValue(binaryOp);
    }
    (*body->body_begin()) = rtnStmt;
    if (!ci.getASTConsumer().HandleTopLevelDecl(DeclGroupRef(add1))) [[unlikely]]
        throw LLVMException(llvm::make_error<llvm::StringError>("[Error]: Consumer rejected the decl", std::error_code()));
    return add1;
}
