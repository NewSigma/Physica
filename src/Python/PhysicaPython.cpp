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
#include <nanobind/stl/string.h>
#include "Physica/Core/Exception/LLVMException.h"
#include "Physica/Python/PhysicaPython.h"
#include "Physica/Python/CXXPtr.h"
#include "Physica/Python/CXXObj.h"
#include "Physica/Python/LLVM/ASTCursor.h"
#include "Physica/Python/LLVM/LLVM.h"

using namespace Physica;
namespace {
    constinit PhysicaPython* instance = nullptr;
}

PhysicaPython::PhysicaPython(std::filesystem::path root_) : clang(std::move(root_), llvm) {
    using namespace llvm::orc;

    const auto& target = clang.getTarget();
    JITTargetMachineBuilder tmBuilder(target.getTriple());
    tmBuilder.addFeatures(target.getTargetOpts().Features);
    LLJITBuilder jitBuilder{};
    jitBuilder.setJITTargetMachineBuilder(std::move(tmBuilder));
    jit = check(jitBuilder.create());

    strTypeMap["NoneType"] = CXXType(ffi_type_void);
    strTypeMap["float"] = CXXType(ffi_type_float);
    strTypeMap["double"] = CXXType(ffi_type_double);
}

void PhysicaPython::compile(const char* moduleName) {
    auto& ptu = clang.compile(moduleName);
    auto& jit = getJIT();
    check(jit.addIRModule(llvm::orc::ThreadSafeModule(std::move(ptu.unitModule), llvm.getThreadSafeContext())));
}

PhysicaPython& PhysicaPython::getInstance() noexcept {
    return *instance;
}

void pymain(nanobind::module_& m);

NB_MODULE(PhysicaPython, m) {
    m.doc() = "Backend of Physica Python interface";
    std::ignore = nanobind::exception<LLVMException>(m, "LLVMException");

    nanobind::class_<CXXPtr>(m, "CXXPtr")
            .def("__repr__", &CXXPtr::toString);

    nanobind::class_<CXXObj>(m, "CXXObj")
            .def(nanobind::init<CXXPtr, nanobind::args>())
            .def("__del__", [](CXXObj& obj) { obj.~CXXObj(); })
            .def("construct", &CXXObj::construct)
            .def("call", &CXXObj::call, nanobind::rv_policy::move);

    nanobind::class_<ASTCursor>(m, "ASTCursor")
            .def(nanobind::new_([]() {
                return ASTCursor(PhysicaPython::getInstance().getClang().getASTContext());
            }))
            .def("__repr__", &ASTCursor::toString)
            .def("pop", &ASTCursor::pop)
            .def("ls", &ASTCursor::ls)
            .def("cd", &ASTCursor::cd)
            .def("dump", &ASTCursor::dump)
            .def("reset", &ASTCursor::reset)
            .def("size", &ASTCursor::size);

    m.def("init", [](const char* root) {
        if (instance == nullptr)
            instance = new PhysicaPython(root);
    });

    m.def("include", [](const char* includeName) -> CXXPtr {
        return (void*)PhysicaPython::getInstance().getClang().include(includeName);
    }, nanobind::rv_policy::move);

    m.def("compile", [](const char* moduleName) {
        PhysicaPython::getInstance().compile(moduleName);
    });

    pymain(m);
}
