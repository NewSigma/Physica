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
#include <Physica/Python/PhysicaPython.h>
#include <Physica/Python/CXXPtr.h>
#include <Physica/Python/CXXObj.h>
#include <Physica/Python/LLVM/LLVM.h>
#include <Physica/Python/Exception/LLVMException.h>

namespace Physica::Python {
    PhysicaPython::PhysicaPython() {
        exec = Executor(clang.getCI().getTarget());

        strTypeMap["NoneType"] = CXXType(ffi_type_void);
        strTypeMap["float"] = CXXType(ffi_type_float);
        strTypeMap["double"] = CXXType(ffi_type_double);
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

    py::class_<CXXPtr>(m, "CXXPtr")
        .def("__repr__", &CXXPtr::toString);

    py::class_<CXXObj>(m, "CXXObj")
        .def(py::init([](CXXPtr p) {
            using namespace clang;
            assert(llvm::isa<CXXRecordDecl>(*p.cast<int>()) && "[Error]: This is not a plain class");
            return CXXObj::create(p.cast<CXXRecordDecl>());
        }))
        .def("__del__", [](CXXObj& obj) { obj.~CXXObj(); })
        .def("call", &CXXObj::call, py::return_value_policy::move);

    m.def("include", [](const char* includeName) -> CXXPtr {
        return (void*)PhysicaPython::getInstance().getClang().include(includeName);
    }, py::return_value_policy::move);

    m.def("compile", [](const char* moduleName) {
        auto& pp = PhysicaPython::getInstance();
        auto* unit = pp.getClang().compile(moduleName);
        if (unit == nullptr)
            return false;
        auto& jit = pp.getJIT();
        llvmCheck(jit.addIRModule(llvm::orc::ThreadSafeModule(std::move(unit->unitModule), LLVM::getInstance().getThreadSafeContext())));
        return true;
    });
}
