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
#include <pybind11/pybind11.h>
#include "Physica/Python/PhysicaPython.h"
#include "Physica/Python/LLVM/LLVM.h"
#include "Physica/Python/LLVM/Clang.h"

using namespace Physica::Python;
namespace py = pybind11;

PYBIND11_MODULE(PhysicaPython, m) {
    m.doc() = "PhysicaPython is a python interface to Physica";

    py::enum_<ScalarOption>(m, "ScalarOption", py::arithmetic())
        .value("Float", ScalarOption::Float)
        .value("Double", ScalarOption::Double)
        .value("MultiPrecision", ScalarOption::MultiPrecision)
        .export_values();

    Class<Scalar<Float>>::pybind(m);
    Class<Scalar<Double>>::pybind(m);
    Class<Physica::Utils::Array<Scalar<Double>>>::pybind(m);
    Class<Vector<Scalar<Double>>>::pybind(m);
    m.def("runKernel", []() -> py::bytes {
        return py::bytes();
    });
    m.def("initLLVM", []() { [[maybe_unused]] auto& llvm = LLVM::getInstance(); });
}
