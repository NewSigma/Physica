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
#include "clang/AST/GlobalDecl.h"
#include "Physica/Python/PhysicaPython.h"
#include "Physica/Python/CXXObj.h"
#include "Physica/Python/LLVM/LLVM.h"
#include "Physica/Python/Exception/LLVMException.h"

namespace Physica::Python {
    PhysicaPython::PhysicaPython() {
        exec = Executor(clang.getCI().getTarget());
    }

    typename PhysicaPython::ExecutorAddr PhysicaPython::emit(const clang::FunctionDecl& func) {
        auto& jit = getJIT();
        auto& ptu = clang.makePTU("Del");
        llvmCheck(jit.addIRModule(llvm::orc::ThreadSafeModule(std::move(ptu.unitModule), LLVM::getInstance().getThreadSafeContext())));
        auto decl = clang::GlobalDecl(&func);
        const std::string symbol = clang.getCodeGen().GetMangledName(std::move(decl)).str();
        return llvmCheck(jit.lookup(symbol));
    }

    PhysicaPython& PhysicaPython::getInstance() noexcept {
        static PhysicaPython instance{};
        return instance;
    }
}

PYBIND11_MODULE(PhysicaPython, m) {
    using namespace Physica::Python;

    m.doc() = "PhysicaPython is a python interface to Physica";
    py::register_exception<LLVMException>(m, "LLVMException");

    py::class_<CXXObj>(m, "CXXObj")
        .def("__del__", [](CXXObj& obj) { obj.~CXXObj(); })
        .def("call", &CXXObj::call, py::return_value_policy::move);

    m.def("include", [](const char* includeName) {
        const void* pClass = PhysicaPython::getInstance().getClang().include(includeName);
        return pClass;
    });
}
